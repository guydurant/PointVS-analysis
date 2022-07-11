from rdkit import Chem
from rdkit.Chem.rdMolTransforms import ComputeCentroid
from vina import Vina
from tqdm import tqdm
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from meeko import PDBQTMolecule
import os 


pdbs = []
with open('pdbbind_2020_general_cleaned_ids.txt', 'r') as f:
    for line in f.readlines():
        pdbs.append(line[:4])
casf_pdbs = []
with open('casf_2016_ids.txt', 'r') as f:
    for line in f.readlines():
        casf_pdbs.append(line[:4])
#missed_pdbs = []
#with open('pdbbind_2020_general_cleaned_ids_w_pymol.txt', 'r') as f:
 #   for line in f.readlines():
  #      if line[:4] not in pdbs:    
   #         pdbs.append(line[:4])

#print(len(missed_pdbs)

v = Vina(sf_name='vina', verbosity=0)
for p in tqdm(casf_pdbs):
    try:
        sdf = f'pdbbind_2020_general/{p}/{p}_ligand.sdf'
        lig = next(Chem.SDMolSupplier(str(sdf), removeHs=False))

    # but we'll try the .mol2 if RDKit can't parse the .sdf
        if lig is None:
            mol2 = f'pdbbind_2020_general/{p}/{p}_ligand.mol2'
            lig = Chem.MolFromMol2File(str(mol2), removeHs=False)
        lig_H = Chem.AddHs(lig)
        AllChem.EmbedMolecule(lig_H) # Regenerate 3D coordinates??

    # prepare ligand pdbqt
        meeko_prep = MoleculePreparation(hydrate=True)
        meeko_prep.prepare(lig_H)
        lig_pdbqt = meeko_prep.write_pdbqt_string()

    

    #RDKit can't read pdbqt format so I use SDF for centroid calculation
        centroid = ComputeCentroid(lig.GetConformer())

        v.set_receptor(f'pdbbind_2020_general/{p}/{p}_protein_cleaned.pdbqt')
        v.set_ligand_from_string(lig_pdbqt)

        v.compute_vina_maps(center=[centroid.x, centroid.y, centroid.z], box_size=[25, 25, 25])

    # Dock the ligand
        v.dock(exhaustiveness=12, n_poses=20)
        output_pdbqt = v.poses(n_poses=20)

    # convert to SDF and write
        pmol = PDBQTMolecule(output_pdbqt)
        f = Chem.SDWriter(f'redocked_pdbbind_2020_general/{p}_redocking-results.sdf')
        for pose in pmol:
            output_rdmol = pmol.export_rdkit_mol()
            output_rdmol_w_bond_order = AllChem.AssignBondOrdersFromTemplate(lig_H, output_rdmol)
            f.write(output_rdmol_w_bond_order)
        f.close()
    except:
        print(f'{p} could not be docked')
