#!/usr/bin/env python3

import sys
import csv
from bs4 import BeautifulSoup
import os

db_file_path = os.path.realpath(os.path.dirname(__file__))+'/metanetx'

reaction_xml = sys.argv[1]
outputsmiles = sys.argv[2]

with open(reaction_xml, 'r') as f:
   reaction_data = f.read()

reaction_parser = BeautifulSoup(reaction_data, "xml")

################## Read MetaNetX ##################

compref = {}
with open(db_file_path+'/chem_xref.tsv', mode='r') as tsv:
   tsv_reader = csv.reader(filter(lambda row: row[0]!='#', tsv), delimiter='\t')
   for splits in tsv_reader:
      compref[splits[0]] = splits[1]

compounds = {}
#ID	name	reference	formula	charge	mass	InChI	InChIKey	SMILES
with open(db_file_path+'/chem_prop.tsv', mode='r') as tsv:
   tsv_reader = csv.reader(filter(lambda row: row[0]!='#', tsv), delimiter='\t')
   for splits in tsv_reader:

      meta_id = splits[0]
      name = splits[1]
      smiles = splits[8]
      
      inchi = splits[7][9:]
      compref[inchi] = meta_id
         
      compounds[meta_id] = (name, smiles)

comp_deprecated = {}
with open(db_file_path+'/chem_depr.tsv', mode='r') as tsv:
   for line in tsv:
      if line.startswith("#"):
         continue
       
      id_old, id_new, version = line.strip().split('\t')
      version = int(version.replace('*', '0').replace('.',''))
     
      if id_old in comp_deprecated and version < comp_deprecated[id_old][1]:
         continue

      if id_old in comp_deprecated:
         comp_deprecated[id_old][0].append(id_new) 
      else:
         comp_deprecated[id_old] = ([id_new], version)

reactions = {}
#ID	mnx_equation	reference	classifs	is_balanced	is_transport
with open(db_file_path+'/reac_prop.tsv', mode='r') as tsv:
   tsv_reader = csv.reader(filter(lambda row: row[0]!='#', tsv), delimiter='\t')
   for splits in tsv_reader:

      meta_id = splits[0]
      mnx_equation = splits[1]
      ecs = [ec for ec in splits[3].split(";") if not ec == ''] 

      reactions[meta_id] = (mnx_equation, ecs)

checkref = {}
with open(db_file_path+'/reac_xref.tsv', mode='r') as tsv:
   tsv_reader = csv.reader(filter(lambda row: row[0]!='#', tsv), delimiter='\t')
   for splits in tsv_reader:

      if ":" not in splits[0]:
         continue

      source, ref_id = splits[0].split(":")
      meta_id = splits[1]

      if source == "bigg.reaction":
         checkref[ref_id] = meta_id

################## Parse Compounds from  XML ##################

xmlcompound_to_potential_metanetxids = {}
for xml_species in reaction_parser.find_all('species'):
   
   xm_id = xml_species.get("id")
   identifiers = []
   identifiers.append("bigg.metabolite:"+xml_species.get("id")[2:-2])
   
   for annotation in xml_species.find_all('rdf:li'):
       resource = annotation.get('rdf:resource')
       resource_split = resource.split("/")[3:]
       if len(resource_split) > 1:
          db, ident = resource_split
          
          if db == "metanetx.chemical":
              identifiers.append(ident)
          elif db == "inchikey":
              identifiers.append(ident)
          elif db == "kegg.compound":
              identifiers.append("kegg.compound:"+ident)
          elif db == "biocyc":
              identifiers.append("metacyc.compound:"+ident.split(":")[1])
          elif db == "seed.compound":
              identifiers.append("seed.compound:"+ident)
       else:
          ident = resource_split[0]
          if ident.startswith("inchikey"):
              identifiers.append(ident.split(":")[1])
          elif ident.startswith("kegg.compound"):
              identifiers.append(ident)
          elif ident.startswith("bigg.metabolite"):
              identifiers.append(ident)
          elif ident.startswith("biocyc:META"):
              identifiers.append("metacyc.compound:"+ident.split(":")[2])

   potential_meta_ids = set()
   for ident in identifiers:
      if ident in compref:
          potential_meta_ids.add(compref[ident])

   while True:
      changed = False
      pmid_it = list(potential_meta_ids)
      for pmid in pmid_it:
         if pmid in comp_deprecated:
            potential_meta_ids.remove(pmid)
            new_list, _ = comp_deprecated[pmid]
            for nlid in new_list:
               potential_meta_ids.add(nlid)
            changed = True
      if not changed:
         break


   xmlcompound_to_potential_metanetxids[xm_id] = {}
   for pmid in potential_meta_ids:
      xmlcompound_to_potential_metanetxids[xm_id][pmid] = 0


################## Find best Entry for each Compound ##################

for xml_react in reaction_parser.find_all('reaction'):

     if xml_react.find('listOfProducts') == None:
        continue # skip sinks

     bigg_id = xml_react.get('id')
     
     if bigg_id in checkref:
        meta_reaction_id = checkref[bigg_id]
        
        reaction_string = reactions[meta_reaction_id][0]
        if reaction_string == " = ":
           continue
        
        meta_ids = []
        reaction_string = reaction_string.replace('=', "+") 
        for rs in reaction_string.split('+'):
           _, comp_id_str = rs.strip().split()
           comp_id = comp_id_str.strip().split('@')[0]
           meta_ids.append(comp_id)
        
        for educt in xml_react.find('listOfReactants').find_all('speciesReference'):
           educt_id = educt.get('species')
           for mid in meta_ids:
              if mid in xmlcompound_to_potential_metanetxids[educt_id]:
                 xmlcompound_to_potential_metanetxids[educt_id][mid] += 1
      
        for product in xml_react.find('listOfProducts').find_all('speciesReference'):
           product_id = product.get('species')
           for mid in meta_ids:
              if mid in xmlcompound_to_potential_metanetxids[product_id]:
                 xmlcompound_to_potential_metanetxids[product_id][mid] += 1       


xmlcompound_to_metanetxid = {}
for c in xmlcompound_to_potential_metanetxids:
   best_meta_id = max(xmlcompound_to_potential_metanetxids[c], key= lambda x: xmlcompound_to_potential_metanetxids[c][x])
   xmlcompound_to_metanetxid[c] = best_meta_id

################## Create SMILES from  XML ##################

def get_compound_info(compound_id):

   name, smiles = compounds[compound_id]
   if smiles.count('.') > 0:
      newname = name + " Subcomponent 1"
      for i in range(smiles.count('.')):
         newname += " + " + name + " Subcomponent "+str(i+2)
      name = newname

   return name, smiles

def parse_and_print_xml_reaction(educts, products):

   name_formula = ""
   smiles_formula = ""

   first = True 
   for xml_id, count in educts:

      if xml_id not in xmlcompound_to_metanetxid:
         return (False, "", "")

      meta_id = xmlcompound_to_metanetxid[xml_id]
      name, smiles = get_compound_info(meta_id)

      if not smiles:
         return (False, "", "")

      for i in range(count):
         if not first:
            name_formula += " + "  
         name_formula += name
      
         if not first:
            smiles_formula += "."
            
         first = False
         smiles_formula += smiles
      
   name_formula += " = "   
   smiles_formula += ">>"
   
   first = True 
   for xml_id, count in products:

      if xml_id not in xmlcompound_to_metanetxid:
         return (False, "", "")

      meta_id = xmlcompound_to_metanetxid[xml_id]
      name, smiles = get_compound_info(meta_id)

      if not smiles:
         return (False, "", "")

      for i in range(count):
         if not first:
            name_formula += " + "  
         name_formula += name
      
         if not first:
            smiles_formula += "."
            
         first = False
         smiles_formula += smiles
      
   return (True, name_formula, smiles_formula)


################## Parse Reactions and Compare ##################

with open(outputsmiles, 'w') as omf:
 for xml_react in reaction_parser.find_all('reaction'):

   if xml_react.find('listOfProducts') == None:
      continue # skip sinks
      
   reversible = xml_react.get('reversible').lower() == "true"
   bigg_id = xml_react.get('id')
   meta_reaction_id = "-"
   if bigg_id in checkref:
      meta_reaction_id = checkref[bigg_id]

   educts = []
   for educt in xml_react.find('listOfReactants').find_all('speciesReference'):
      educts.append( (educt.get('species'), int(float(educt.get('stoichiometry')))) )
      
   products = []
   for product in xml_react.find('listOfProducts').find_all('speciesReference'):
      products.append( (product.get('species'), int(float(product.get('stoichiometry')))) )
      

   valid, name_formula, smiles_formula = parse_and_print_xml_reaction(educts, products)
   if valid:
      ec_string = "-"
      if meta_reaction_id != "-":
         ec_string = ";".join(reactions[meta_reaction_id][1])
      print(file=omf)
      print("Bigg ID:", bigg_id, "MetaNetXId:", meta_reaction_id, "Reversible:", reversible, file=omf)
      print("ECs:", ec_string, file=omf)
      print(name_formula, file=omf)   
      print(smiles_formula, file=omf)





   

