#!/bin/bash   
#SBATCH -J attribution                     # Job name
#SBATCH -A durant                         # Project Account                                   
#S BATCH --time=00:30:00                 # Walltime                                      
#SBATCH --mem-per-cpu=4096             # memory/cpu (in MB) ### commented out              
#SBATCH --ntasks=1                      # 1 tasks                                               
#SBATCH --cpus-per-task=2               # number of cores per task                          
#SBATCH --nodes=1 # number of nodes
#SBATCH --gres=gpu:1
#SBATCH --overcommit
#SBATCH --gpus-per-node=1
#S BATCH --exclusive                     # node should not be shared with other jobs, only use this if you intend the node to be usable only by you as this will block other users from submitting jobs to the same node                
#SBATCH --chdir=/data/localhost/not-backed-up/durant/PointVS # From where you want the job to be run
#SBATCH --mail-user=guy.durant@linacre.ox.ac.uk  # set email address                           
#S BATCH --mail-type=ALL                 # Spam us with everything, caution
#SBATCH --mail-type=begin               # Instead only email when job begins...
#SBATCH --mail-type=end                 # ... and ends
#SBATCH --clusters=swan
#SBATCH --partition=high-opig-gpu    # Select a specific partition rather than default 
#S BATCH -w naga03.cpu.stats.ox.ac.uk # Provide a specific node/nodelist rather than the standard nodelist associated with the partition (useful if you have a data setup on one specific node)
#SBATCH --output=/data/localhost/not-backed-up/durant/slurm_%j.out  # Writes standard output to this file. %j is jobnumber                             
#SBATCH --error=/data/localhost/not-backed-up/durant/slurm_%j.err   # Writes error messages to this file. %j is jobnumber

source ~/.bashrc
source ~/miniconda3/etc/profile.d/conda.sh
conda activate pointvs
echo The batch script has been started

python point_vs/attribution/attribution.py atom_masking crystal_regression/PointBAP/crystal_regression_both_30/checkpoints/ckpt_epoch_9.pt 4abg_attribution_unbiased --output_name 4abg_attribution_unbiased_crystal --input_receptor_file ../pdbbind_2020_general/4abg/4abg_protein_cleaned.pdb --input_ligand_file ../pdbbind_2020_general/4abg/4abg_ligand.sdf

python point_vs/attribution/attribution.py atom_masking crystal_regression/PointBAP/crystal_regression_both_30/checkpoints/ckpt_epoch_9.pt 4abg_attribution_unbiased --output_name 4abg_attribution_unbiased_redocked --input_receptor_file ../pdbbind_2020_general/4abg/4abg_protein_cleaned.pdb --input_ligand_file ../pdbbind_2020_general/4abg/4abg_redocking_best_pose.sdf

python point_vs/attribution/attribution.py atom_masking crystal_regression/PointBAP/crystal_regression_both_30/checkpoints/ckpt_epoch_9.pt 4abg_attribution_unbiased --output_name 4abg_attribution_unbiased_crossdocked_most --input_receptor_file ../pdbbind_2020_general/1k1i/1k1i_protein_cleaned.pdb --input_ligand_file ../pdbbind_2020_general/4abg/4abg_crossdocking_most_structure_similiar_best_pose.sdf

python point_vs/attribution/attribution.py atom_masking crystal_regression/PointBAP/crystal_regression_both_30/checkpoints/ckpt_epoch_9.pt 4abg_attribution_unbiased --output_name 4abg_attribution_unbiased_crosdocked_least --input_receptor_file ../pdbbind_2020_general/1uto/1uto_protein_cleaned.pdb --input_ligand_file ../pdbbind_2020_general/4abg/4abg_crossdocking_least_structure_similiar_best_pose.sdf

python point_vs/attribution/attribution.py atom_masking crystal_regression/PointBAP/crystal_regression/checkpoints/ckpt_epoch_9.pt 4abg_attribution --output_name 4abg_attribution_biased_crystal --input_receptor_file ../pdbbind_2020_general/4abg/4abg_protein_cleaned.pdb --input_ligand_file ../pdbbind_2020_general/4abg/4abg_ligand.sdf

python point_vs/attribution/attribution.py atom_masking crystal_regression/PointBAP/crystal_regression/checkpoints/ckpt_epoch_9.pt 4abg_attribution --output_name 4abg_attribution_biased_redocked --input_receptor_file ../pdbbind_2020_general/4abg/4abg_protein_cleaned.pdb --input_ligand_file ../pdbbind_2020_general/4abg/4abg_redocking_best_pose.sdf

python point_vs/attribution/attribution.py atom_masking crystal_regression/PointBAP/crystal_regression/checkpoints/ckpt_epoch_9.pt 4abg_attribution --output_name 4abg_attribution_biased_crossdocked_most --input_receptor_file ../pdbbind_2020_general/1k1i/1k1i_protein_cleaned.pdb --input_ligand_file ../pdbbind_2020_general/4abg/4abg_crossdocking_most_structure_similiar_best_pose.sdf

python point_vs/attribution/attribution.py atom_masking crystal_regression/PointBAP/crystal_regression/checkpoints/ckpt_epoch_9.pt 4abg_attribution --output_name 4abg_attribution_biased_crossdocked_least  --input_receptor_file ../pdbbind_2020_general/1uto/1uto_protein_cleaned.pdb --input_ligand_file ../pdbbind_2020_general/4abg/4abg_crossdocking_least_structure_similiar_best_pose.sdf













echo The script has finised running
