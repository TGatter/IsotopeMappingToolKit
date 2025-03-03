#!/usr/bin/env python3

import os
import sys
import enum

import logging

from custom_pysmiles import read_smiles
from custom_pysmiles.smiles_helper import (add_explicit_hydrogens, remove_explicit_hydrogens)

import networkx as nx
from networkx.algorithms import isomorphism as nxisomorphism

from rdkit import Chem

from pyvis.network import Network

#logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s', level=logging.INFO)
logging.basicConfig(format='%(levelname)s:\t%(message)s', level=logging.DEBUG)

NO_MAP_DEFAULT_KEY = -1

mappedsmiles = sys.argv[1]
outputgml = sys.argv[2] 

molecule_graph_path = sys.argv[3]

map_hydrogens = False
if len(sys.argv) == 5:
    map_hydrogens = str(sys.argv[4])
    
# generate Folder with all molecules as single graphs
newpath = 'Molecule_Graphs' 
if not os.path.exists(molecule_graph_path):
    os.makedirs(molecule_graph_path)

# reading in the list of highly concentrated molecules 
with open ( 'metanetx/list_highlyConcMol.txt' , 'r') as highmol_file:
    highmol_list = [s.strip() for s in highmol_file.readlines() ]


## ========== DEF BLOCK

@enum.unique
class TransitionType(str, enum.Enum):
    """Possible SMILES token types"""
    NO_TRANSITION = "ChemicalBond"
    SYMMETRY = "Symmetry"
    REACTION = "Reaction"
    HYDROGEN_GROUP = "HydrogenGroup"
    HYDROGEN_REACTION = "HydrogenReaction"
    HYDROGEN_FREE = "HydrogenFreedReaction"

@enum.unique
class DirectionType(enum.IntEnum):
    """Possible SMILES token types"""
    UNDIRECTED = 0
    BIDIRECTED = 1
    DIRECTED = 2

def hasPath(G, s, t, edge_type):
    if s == t:
        return True

    visited = set()
    visited.add(s)
    stack = [(s, iter(G[s]))]
    while stack:
        parent, children = stack[-1]

        for child in children:
            if child not in visited and G.edges[parent, child]['transition'] == edge_type:
                if child == t:
                    return True
                visited.add(child)
                stack.append((child, iter(G[child])))
        stack.pop()
    return False

def findHydrogenGroups(mol, mapped_hydrogens, mapped_atoms):

    inv_mapped_atom = {v: k for k, v in mapped_atoms.items()}    

    for n in mol.nodes():
        if mol.nodes[n].get('element', '') == 'H':
            class_id = NO_MAP_DEFAULT_KEY
            for neighbor in mol[n]:
                if neighbor in inv_mapped_atom:
                    class_id = inv_mapped_atom[neighbor]
            mapped_hydrogens.setdefault(class_id, []).append(n)

def findIsomorphATNStructure(ATN, mol):

    em = nxisomorphism.categorical_edge_match(['order'],[0])
    nm = nxisomorphism.categorical_node_match(['element', 'isotope', 'hcount', 'charge'],['', 0, 0, 0])
    GM = nxisomorphism.GraphMatcher(ATN, mol, node_match=nm, edge_match=em)   
    result_iso = GM.subgraph_is_monomorphic()
    mapping_ATN_to_mol = GM.mapping   

    return (result_iso, mapping_ATN_to_mol)

def addAutomorphisms(mol, limit_to_orbits=True):
    em = nxisomorphism.categorical_edge_match(['order'],[0])
    nm = nxisomorphism.categorical_node_match(['element', 'isotope', 'hcount', 'charge'],['', 0, 0, 0])
    GM = nxisomorphism.GraphMatcher(mol, mol, node_match=nm, edge_match=em)   

    permutation_lists = list(GM.isomorphisms_iter())

    if limit_to_orbits:
        blockset = set()
        for node in mol.nodes():
            if node in blockset:
                continue
            for isomorphism in permutation_lists:
                if node != isomorphism[node]:
                   mol.add_edge(node, isomorphism[node])
                   mol.edges[node, isomorphism[node]]['transition'] = TransitionType.SYMMETRY
                   mol.edges[node, isomorphism[node]]['directed'] = DirectionType.BIDIRECTED
                   blockset.add(isomorphism[node])
    else:
        for isomorphism in permutation_lists:
            for i in isomorphism:
                if i != isomorphism[i] and not mol.has_edge(i, isomorphism[i]):
                    mol.add_edge(i, isomorphism[i])
                    mol.edges[i, isomorphism[i]]['transition'] = TransitionType.SYMMETRY
                    mol.edges[i, isomorphism[i]]['directed'] = DirectionType.BIDIRECTED
       
def parseXDuct(name, smiles, compound_to_subgraph, compoundId_to_compound, ATN, mapped_atoms, mapped_hydrogens = {}, explicit_hydrogens=False):

    logging.debug("Parse " + name + " : " + smiles)
    
    fixed_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), allHsExplicit=True, canonical=False)
    mol = read_smiles( fixed_smiles ) # read smile

    # RXNMApper does not create valid H mappings, remove them all
    for n in mol.nodes():
        if mol.nodes[n].get('element', '') == 'H' and 'class' in mol.nodes[n]:
            del mol.nodes[n]['class']
    remove_explicit_hydrogens(mol)
    
    if name in compound_to_subgraph:

        logging.debug("Use Existing " + name)

        if explicit_hydrogens:
            add_explicit_hydrogens(mol)

        # we already added this compound to the network
        compound_subgraph = ATN.subgraph(compound_to_subgraph[name])
        has_isomorph_subgraph, mapping_ATN_to_mol = findIsomorphATNStructure(compound_subgraph, mol)
        if has_isomorph_subgraph:
            for atn_node in mapping_ATN_to_mol:
                if 'class' in mol.nodes[mapping_ATN_to_mol[atn_node]]:
                    mapped_atoms[ mol.nodes[mapping_ATN_to_mol[atn_node]]['class'] ] = atn_node
        else:
            logging.error("Compound Naming Error: " + name + " : " + smiles)

        if explicit_hydrogens:
            logging.debug("Find Hydrogen Partners")
            findHydrogenGroups(compound_subgraph, mapped_hydrogens, mapped_atoms)

    else:

        logging.debug("Add New Compound "+ name)

        # this is definitely new
        nextCId = len(compound_to_subgraph)
        compoundId_to_compound[nextCId] = (name, fixed_smiles)

        rename = {node : str(nextCId)+'_'+str(node) for node in mol.nodes()} # rename all nodes so that we cannot have collisions in the ATN
        nx.relabel_nodes(mol, rename, copy=False)

        addAutomorphisms(mol)

        if explicit_hydrogens:
            add_explicit_hydrogens(mol, str(nextCId)+"_")

        for n in mol.nodes():
            mol.nodes[n]['compound_id'] = nextCId
            mol.nodes[n]['compound_name'] = name

        for e in mol.edges():
            mol.edges[e]['compound_id'] = nextCId
            mol.edges[e]['compound_name'] = name
            if 'transition' not in mol.edges[e]:
                mol.edges[e]['transition'] = TransitionType.NO_TRANSITION
                mol.edges[e]['directed'] = DirectionType.UNDIRECTED

        ATN.add_nodes_from(mol.nodes(data=True))
        ATN.add_edges_from(mol.edges(data=True))

        for mol_node, data in mol.nodes(data=True):
            if 'class' in data:
                mapped_atoms[data['class']] = mol_node
                del mol.nodes[mol_node]['class']
            compound_to_subgraph.setdefault(name, []).append(mol_node)

        if explicit_hydrogens:
            logging.debug("Find Hydrogen Partners")
            findHydrogenGroups(mol, mapped_hydrogens, mapped_atoms)
            
        
# ======== MAIN

# first an intermediary with bonding and transition edges that is trimmed afterwards
ATN=nx.Graph()
compound_to_subgraph = {}
compoundId_to_compound = {}

reactions = []

# load all SMILES from txt files
with open( mappedsmiles , 'r') as smiles_file:
  while True:    
    meta_line = smiles_file.readline()
    name_line = smiles_file.readline()
    mapped_smiles_line = smiles_file.readline()
    if not name_line or not  mapped_smiles_line:
       break;

    logging.info("Next Reaction ==============")
    logging.info(name_line.strip())
    logging.info(mapped_smiles_line.strip())

    smiles_str = mapped_smiles_line.strip().replace("@", '').replace("/", '')
    names_str = name_line.strip()
    reversable = meta_line.strip().split()[6] == "True"

    #filename = ''.join(letter for letter in name_line if letter.isalnum())

    reactions.append( (names_str, smiles_str) )

    # create left and right list with the smiles
    smiles_sides = smiles_str.split('>>')
    smiles_left = smiles_sides[0].split('.')
    smiles_right = smiles_sides[1].split('.')
    
    # create left and right list with the metabolite names
    names_sides = names_str.split('=')
    names_left = list(map(str.strip, names_sides[0].split(' + ')))
    names_right = list(map(str.strip, names_sides[1].split(' + ')))

    # we have compartment change reactions that are not interesting for the ATN 
    if set(names_left) == set(names_right):
        continue

    logging.debug("Parse Educts")
    mapped_educt_atoms = {}
    mapped_educt_hydrogens = {}
    for name, smiles in zip(names_left, smiles_left):

        modifier = ''  # rename molecules of highly conzentrations
        if name in highmol_list:
            modifier = '_in'

        parseXDuct(name+modifier, smiles, compound_to_subgraph, compoundId_to_compound, ATN, mapped_educt_atoms, mapped_educt_hydrogens, map_hydrogens)

    logging.debug("Parse Products")
    mapped_product_atoms = {}
    mapped_product_hydrogens = {}
    for name, smiles in zip(names_right, smiles_right):
        
        modifier = ''  # rename molecules of highly conzentrations
        if name in highmol_list:
            modifier = '_out'
    
        parseXDuct(name+modifier, smiles, compound_to_subgraph, compoundId_to_compound, ATN, mapped_product_atoms, mapped_product_hydrogens, map_hydrogens)

    if map_hydrogens:

        logging.debug("Map Hydrogens")

        trans_H_node = "react_"+str(len(reactions)-1)+"_free_H"

        key_set_educt  = set(mapped_educt_hydrogens.keys())
        key_set_product  = set(mapped_product_hydrogens.keys())
        key_set_all = key_set_educt.union(key_set_product)

        for key in key_set_all:

            hydro_count_educts = 0
            if key in key_set_educt:
                hydro_count_educts = len(mapped_educt_hydrogens[key])
            hydro_count_products = 0
            if key in key_set_product:
                hydro_count_products = len(mapped_product_hydrogens[key])

            if key == NO_MAP_DEFAULT_KEY or not hydro_count_educts == hydro_count_products:
                if not ATN.has_node(trans_H_node):
                    ATN.add_node(trans_H_node)
                    ATN.nodes[trans_H_node]['element'] = 'H'
                    
                if hydro_count_educts > hydro_count_products or (key == NO_MAP_DEFAULT_KEY and hydro_count_educts > 0):
                    for n in mapped_educt_hydrogens[key]:
                        ATN.add_edge(n, trans_H_node)
                        ATN.edges[n, trans_H_node]['reaction_id'] = {}
                        ATN.edges[n, trans_H_node]['directed'] = DirectionType.DIRECTED

                        ATN.edges[n, trans_H_node]['transition'] = TransitionType.HYDROGEN_FREE
                        ATN.edges[n, trans_H_node]['moving_atom'] = hydro_count_educts - hydro_count_products
                        ATN.edges[n, trans_H_node]['reaction_id'][(n, trans_H_node)] = str(len(reactions)-1)

                        if reversable:
                            ATN.edges[n, trans_H_node]['directed'] = DirectionType.BIDIRECTED
                            ATN.edges[n, trans_H_node]['reaction_id'][(trans_H_node, n)] = str(len(reactions)-1)

                if hydro_count_educts < hydro_count_products or (key == NO_MAP_DEFAULT_KEY and hydro_count_products > 0):
                    for n in mapped_product_hydrogens[key]:
                        ATN.add_edge(trans_H_node, n)
                        ATN.edges[trans_H_node, n]['reaction_id'] = {}
                        ATN.edges[trans_H_node, n]['directed'] = DirectionType.DIRECTED

                        ATN.edges[trans_H_node, n]['transition'] = TransitionType.HYDROGEN_FREE
                        ATN.edges[trans_H_node, n]['moving_atom'] = hydro_count_products - hydro_count_educts
                        ATN.edges[trans_H_node, n]['reaction_id'][(trans_H_node, n)] = str(len(reactions)-1)

                        if reversable:
                            ATN.edges[trans_H_node, n]['directed'] = DirectionType.BIDIRECTED
                            ATN.edges[trans_H_node, n]['reaction_id'][(n, trans_H_node)] = str(len(reactions)-1)

            # we always need to add in the real H transitions
            if key == NO_MAP_DEFAULT_KEY or hydro_count_educts == 0 or hydro_count_products == 0:
                continue
 
            rep_atom_educt = mapped_educt_hydrogens[key][0]
            rep_atom_product = mapped_product_hydrogens[key][0]
            if not ATN.has_edge(rep_atom_educt,rep_atom_product):
                ATN.add_edge(rep_atom_educt,rep_atom_product)
                ATN.edges[rep_atom_educt,rep_atom_product]['transition'] = TransitionType.HYDROGEN_REACTION
                ATN.edges[rep_atom_educt,rep_atom_product]['reaction_id'] = {}

            if (rep_atom_educt,rep_atom_product) in ATN.edges[rep_atom_educt,rep_atom_product]['reaction_id']:
                ATN.edges[rep_atom_educt,rep_atom_product]['reaction_id'][(rep_atom_educt,rep_atom_product)] += ","+str(len(reactions)-1)
            else:
                ATN.edges[rep_atom_educt,rep_atom_product]['reaction_id'][(rep_atom_educt,rep_atom_product)] = str(len(reactions)-1)

            if reversable:
                if (rep_atom_product, rep_atom_educt) in ATN.edges[rep_atom_educt,rep_atom_product]['reaction_id']:
                    ATN.edges[rep_atom_educt,rep_atom_product]['reaction_id'][(rep_atom_product, rep_atom_educt)] += ","+str(len(reactions)-1)
                else:
                    ATN.edges[rep_atom_educt,rep_atom_product]['reaction_id'][(rep_atom_product, rep_atom_educt)] = str(len(reactions)-1)          

            for atom in mapped_educt_hydrogens[key][1:]: # all in group need to be connected
                if not hasPath(ATN, rep_atom_educt, atom, TransitionType.HYDROGEN_GROUP):
                    ATN.add_edge(rep_atom_educt, atom)
                    ATN.edges[rep_atom_educt, atom]['transition'] = TransitionType.HYDROGEN_GROUP

            for atom in mapped_product_hydrogens[key][1:]: # all in group need to be connected
                if not hasPath(ATN, rep_atom_product, atom, TransitionType.HYDROGEN_GROUP):
                    ATN.add_edge(rep_atom_product, atom)
                    ATN.edges[rep_atom_product, atom]['transition'] = TransitionType.HYDROGEN_GROUP
  
    for c in mapped_educt_atoms:

        n1 = mapped_educt_atoms[c]
        n2 = mapped_product_atoms[c]

        if ATN.nodes[n1]['compound_id'] == ATN.nodes[n2]['compound_id']:
             continue # skip all self maps

        if not ATN.has_edge(n1,n2):
            ATN.add_edge(n1, n2)
            ATN.edges[n1, n2]['reaction_id'] = {}
            ATN.edges[n1, n2]['transition'] = TransitionType.REACTION
            
        if (n1, n2) in ATN.edges[n1, n2]['reaction_id']:
            ATN.edges[n1, n2]['reaction_id'][(n1, n2)] += ","+str(len(reactions)-1)
        else:
            ATN.edges[n1, n2]['reaction_id'][(n1, n2)] = str(len(reactions)-1)
        
        if reversable and not (ATN.nodes[n2]['compound_name'].endswith('_in') or ATN.nodes[n2]['compound_name'].endswith('_out') or ATN.nodes[n1]['compound_name'].endswith('_in') or ATN.nodes[n1]['compound_name'].endswith('_out')):
            if (n2, n1) in ATN.edges[n1, n2]['reaction_id']:
                ATN.edges[n1, n2]['reaction_id'][(n2, n1)] += ","+str(len(reactions)-1)
            else:
                ATN.edges[n1, n2]['reaction_id'][(n2, n1)] = str(len(reactions)-1)            
 
 
for e in ATN.edges():
    if ATN.edges[e]['transition'] == TransitionType.REACTION or ATN.edges[e]['transition'] == TransitionType.HYDROGEN_REACTION:
        if len(ATN.edges[e]['reaction_id']) == 1:
            ATN.edges[e]['directed'] = DirectionType.DIRECTED
        else:
            ATN.edges[e]['directed'] = DirectionType.BIDIRECTED
   

def add_meta(e, ATN, DATN):
    DATN.edges[e]['transition'] = ATN.edges[e]['transition'].value
    if 'oder' in ATN.edges[e]:
        DATN.edges[e]['order'] = ATN.edges[e]['order']
    if 'compound_id' in ATN.edges[e]:
        DATN.edges[e]['compound_id'] = ATN.edges[e]['compound_id']
    if 'compound_name' in ATN.edges[e]:
        DATN.edges[e]['compound_name'] = ATN.edges[e]['compound_name']
    if 'moving_atom' in ATN.edges[e]:
        DATN.edges[e]['moving_atom'] = ATN.edges[e]['moving_atom']
    if 'reaction_id' in ATN.edges[e] and e in ATN.edges[e]['reaction_id']:
        DATN.edges[e]['reaction_ids'] = ATN.edges[e]['reaction_id'][e]

def bfs_ready_transform(ATN, filter_bonds=True):
    DATN = nx.DiGraph()
    DATN.add_nodes_from(ATN.nodes(data=True))

    for n in DATN.nodes():
        if 'class' in  DATN.nodes[n]:
            del DATN.nodes[n]['class']
        if 'element' not in DATN.nodes[n]:
            DATN.nodes[n]['element'] = "*"

    for e in ATN.edges(): 
        if filter_bonds and ATN.edges[e]['transition'] == TransitionType.NO_TRANSITION:
            continue

        if ATN.edges[e]['directed'] == DirectionType.DIRECTED:
            s,t = list(ATN.edges[e]['reaction_id'].keys())[0] 
            DATN.add_edge(s,t)
            add_meta((s,t), ATN, DATN)
        elif ATN.edges[e]['directed'] == DirectionType.BIDIRECTED or ATN.edges[e]['directed'] == DirectionType.UNDIRECTED:
            s,t = e
            DATN.add_edge(s,t)
            add_meta((s,t), ATN, DATN)
            DATN.add_edge(t,s)
            add_meta((t,s), ATN, DATN)

    return DATN 

nx.write_gml(bfs_ready_transform(ATN), outputgml)

with open(outputgml+".ckey", 'w') as og:
  for key in sorted(compoundId_to_compound):
      print(key, compoundId_to_compound[key][0], compoundId_to_compound[key][1] ,sep='\t', file=og)

# DRAWING
def draw_single_compound(name, mol):

        newMol = mol.copy()
        for n in mol.nodes():
            if 'element' in mol.nodes[n]:
               newMol.nodes[n]['label'] = mol.nodes[n]['element'] + '(' + n + ')'
               if mol.nodes[n]['element'] == 'C':
                  newMol.nodes[n]['color'] = "black"
               if mol.nodes[n]['element'] == 'O':
                  newMol.nodes[n]['color'] = "red"    
               if mol.nodes[n]['element'] == 'C':
                  newMol.nodes[n]['color'] = "black"
               if mol.nodes[n]['element'] == 'S':
                  newMol.nodes[n]['color'] = "yellow"
            else:
                newMol.nodes[n]['label'] = "*"

        for e in mol.edges():
            del newMol.edges[e]['directed']
            if mol.edges[e]['transition'] == TransitionType.REACTION:
                newMol.edges[e]['label'] = str(mol.edges[e]['order'])
                newMol.edges[e]['color'] = "black"
            elif mol.edges[e]['transition'] == TransitionType.SYMMETRY:
                newMol.edges[e]['color'] = "green"


        net = Network()
        net.from_nx(newMol)
        with open(molecule_graph_path + '/' + name + '.html', "w+") as out:
            out.write(net.generate_html())

for comp in compound_to_subgraph:
    compound_subgraph = ATN.subgraph(compound_to_subgraph[comp])        
    draw_single_compound(comp, compound_subgraph)

draw = nx.DiGraph()
draw.add_nodes_from(ATN.nodes())
draw.add_edges_from(ATN.edges(data=True))

for n in ATN.nodes():
    if 'element' in ATN.nodes[n]:
        draw.nodes[n]['label'] = ATN.nodes[n]['element']
    else:
        draw.nodes[n]['label'] = "*"
    if 'hcount' in ATN.nodes[n]:
        draw.nodes[n]['label'] += str(ATN.nodes[n]['hcount'])+'H'
    draw.nodes[n]['color'] = "black"

for e in ATN.edges():

    del draw.edges[e]['directed']
    if 'reaction_id' in draw.edges[e]:
        del draw.edges[e]['reaction_id']
    
    if ATN.edges[e]['directed'] == DirectionType.DIRECTED:
        s,t = list(ATN.edges[e]['reaction_id'].keys())[0] 
        if s==e[0] and t==e[1]:
            draw.edges[e]['arrows'] = "to"
        else:
            draw.edges[e]['arrows'] = "from"
    elif ATN.edges[e]['directed'] == DirectionType.BIDIRECTED:
        draw.edges[e]['arrows'] = "to, from"    
    elif ATN.edges[e]['directed'] == DirectionType.UNDIRECTED:
        draw.edges[e]['arrows'] = "no"

    if ATN.edges[e]['transition'] == TransitionType.SYMMETRY:
        draw.edges[e]['color'] = "green"
    elif ATN.edges[e]['transition'] == TransitionType.REACTION:
        draw.edges[e]['color'] = "red"
    elif ATN.edges[e]['transition'] == TransitionType.HYDROGEN_GROUP or ATN.edges[e]['transition'] == TransitionType.HYDROGEN_REACTION:
        draw.edges[e]['color'] = "LightSkyBlue"
    elif ATN.edges[e]['transition'] == TransitionType.HYDROGEN_FREE:
        draw.edges[e]['color'] = "DarkBlue"
    else:
        draw.edges[e]['label'] = str(ATN.edges[e]['order'])
        draw.edges[e]['color'] = "black"

nt = Network('1000px', '1000px', directed=True)
nt.show_buttons()
nt.from_nx(draw)
with open(molecule_graph_path + '/full_ATN.html', "w+") as out:
    out.write(nt.generate_html())


