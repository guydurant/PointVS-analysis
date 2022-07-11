from rdkit import Chem
from rdkit.Chem.rdMolTransforms import ComputeCentroid
from vina import Vina
from tqdm import tqdm
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from meeko import PDBQTMolecule
import os
import pandas as pd


casf_pdbs = []
with open('casf_2016_ids.txt') as f:
    for line in f:
        casf_pdbs.append(line[:4])

casf_clusters = dict(pd.read_csv('casf_2016_clusters.csv', index_col='PDB code').T)

def append_value(dict_obj, key, value):
    # Check if key exist in dict or not
    if key in dict_obj:
        # Key exist in dict.
        # Check if type of value of key is list or not
        if not isinstance(dict_obj[key], list):
            # If type is not list then make it list
            dict_obj[key] = [dict_obj[key]]
        # Append the value in list
        dict_obj[key].append(value)
    else:
        # As key is not in dict,
        # so, add key-value pair
        dict_obj[key] = value

cluster_dict = {}

for p in casf_pdbs:
    append_value(cluster_dict, casf_clusters[p.upper()]['Cluster ID'], p)

v = Vina(sf_name='vina', verbosity=0)


def prepare_lig_for_docking(pdb, return_only_centroid=False, return_template=False):
    sdf = f'pdbbind_2020_general/{pdb}/{pdb}_ligand.sdf'
    lig = next(Chem.SDMolSupplier(str(sdf), removeHs=False))

    # but we'll try the .mol2 if RDKit can't parse the .sdf
    if lig is None:
        mol2 = f'pdbbind_2020_general/{pdb}/{pdb}_ligand.mol2'
        lig = Chem.MolFromMol2File(str(mol2), removeHs=False)
    lig_H = Chem.AddHs(lig)
    if return_template:
        return lig_H
    AllChem.EmbedMolecule(lig_H) # Regenerate 3D coordinates?
    if return_only_centroid:
        return ComputeCentroid(lig.GetConformer())
    else:
        # prepare ligand pdbqt
        meeko_prep = MoleculePreparation(hydrate=True)
        meeko_prep.prepare(lig_H)
        lig_pdbqt = meeko_prep.write_pdbqt_string()
        return lig_pdbqt


for i in cluster_dict.keys():
    for c in cluster_dict[i]:
        centroid = prepare_lig_for_docking(c, return_only_centroid=True)
        for d in [d for d in cluster_dict[i] if d != c]:
            #try:

            lig_pdbqt = prepare_lig_for_docking(d)
            v.set_receptor(f'pdbbind_2020_general/{c}/{c}_protein_cleaned.pdbqt')
            v.set_ligand_from_string(lig_pdbqt)

            v.compute_vina_maps(center=[centroid.x, centroid.y, centroid.z], box_size=[25, 25, 25])

            # Dock the ligand
            v.dock(exhaustiveness=12, n_poses=20)
            output_pdbqt = v.poses(n_poses=20)

            # convert to SDF and write
            pmol = PDBQTMolecule(output_pdbqt)
            f = Chem.SDWriter(f'crossdocks_casf_2016/{d}_crossdocking_into_{c}-results.sdf')
            for pose in pmol:
                output_rdmol = pmol.export_rdkit_mol()
                lig_H = prepare_lig_for_docking(d, return_template=True)
                output_rdmol_w_bond_order = AllChem.AssignBondOrdersFromTemplate(lig_H, output_rdmol)
                print('file', output_rdmol)
                f.write(output_rdmol_w_bond_order)
            f.close()
            #except:
               # print(f'Could not dock {d} into {c} receptor')
