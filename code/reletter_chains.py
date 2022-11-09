from Bio.PDB import PDBParser, PDBIO
import json
import os

allele_list_path = 'data/alleles.json'
basepath = 'predictions/class_i/raw'
new_basepath = 'predictions/class_i/renumbered'

def change_colabfold_chains(structure):
    for chain in structure.get_chains():
        if chain.id == 'C':
            chain.id = 'A'
    return structure

def reverse_chain_labels(structure):
    for chain in structure.get_chains():
        if chain.id == 'B':
            chain.id = 'D'
        if chain.id == 'A':
            chain.id = 'C'
    for chain in structure.get_chains():
        if chain.id == 'C':
            chain.id = 'B'
        if chain.id == 'D':
            chain.id = 'A' 
    return structure


def change_esmfold_chains(structure):
    structure = reverse_chain_labels(structure)  
    for chain in structure.get_chains():
        if chain.id == 'A':
            print ('renumbering chain A')
            i = 1
            for residue in chain:
                if ' ' in residue.id[2]:
                    residue.id = (residue.id[0], i, residue.id[2])
                    i += 1
    return structure


def change_omgfold_chains(structure):
    structure = reverse_chain_labels(structure)
    return structure


def create_directories(directory_path):
    try:
        os.makedirs(directory_path)
        print (f'{directory_path} created')
    except FileExistsError:
        print (f'{directory_path} already exists')
    pass
        


def build_filepaths(allele, locus, method_slug):
    if method_slug == 'col':
        filepath = f'{basepath}/colabfold/{locus}/{allele}_col/{allele}_col_unrelaxed_rank_1_model_1.pdb'
        method_name = 'colabfold'
    elif method_slug == 'esm':
        filepath = f'{basepath}/esmfold/{locus}/{allele}_esm/{allele}_esmfold.pdb'
        method_name = 'esmfold'
    elif method_slug == 'omg':
        filepath = f'{basepath}/omegafold/{locus}/{allele}_omg/{allele}_omg.pdb'
        method_name = 'omegafold'

    directory_path = f'{new_basepath}/{method_name}/{locus}/{allele}_{method_slug}'
    create_directories(directory_path)
    renumbered_filepath = f'{directory_path}/{allele}_{method_slug}.pdb'
    return filepath, renumbered_filepath


def standardise_prediction(allele, locus, method):
    filepath, renumbered_filepath = build_filepaths(allele, locus, method)
    parser = PDBParser(PERMISSIVE=0)
    structure = parser.get_structure(allele, filepath)

    print (f'{allele=}')
    print (f'{method=}')

    print ('before relettering')
    for chain in structure.get_chains():
        print (chain.id)
        print (len([residue for residue in chain]))


    if method == 'col':
        structure = change_colabfold_chains(structure)
    elif method == 'esm':
        structure = change_esmfold_chains(structure)
    elif method == 'omg':
        structure = change_omgfold_chains(structure)


    print ('after relettering')
    for chain in structure.get_chains():
        print (chain.id)
        print (len([residue for residue in chain]))

    io=PDBIO()
    io.set_structure(structure)
    io.save(renumbered_filepath)



method = 'omg'
allele = 'hla_a_01_01'
locus = 'hla_a'

with open(allele_list_path) as allele_list_file:
    allele_list = json.load(allele_list_file)



for locus in allele_list:
    alleles = allele_list[locus]
    for allele in alleles:
        print ('--------')
        print (allele)
        for method in ['omg','esm','col']:
            standardise_prediction(allele, locus, method)
        print ('--------')















