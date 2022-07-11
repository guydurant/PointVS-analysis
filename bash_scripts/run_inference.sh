#!/bin/bash   
#SBATCH -J inference                     # Job name
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

#python point_vs/inference.py ../crystal_regression/PointBAP/crystal_regression/checkpoints/ckpt_epoch_8.pt /data/localhost/not-backed-up/durant/casf_2016_pymol_redocked.types /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets  --wandb_project PointBAP --wandb_run crystal_regression_redocked_casf
#python point_vs/inference.py ../crystal_regression/PointBAP/crystal_regression/checkpoints/ckpt_epoch_8.pt /data/localhost/not-backed-up/durant/casf_2016_pymol_crossdocked_most_similiar_structure.types /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets  --wandb_project PointBAP --wandb_run crystal_regression_crossdocked_most_casf
#python point_vs/inference.py ../crystal_regression/PointBAP/crystal_regression/checkpoints/ckpt_epoch_8.pt /data/localhost/not-backed-up/durant/casf_2016_pymol_crossdocked_least_similiar_structure.types /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets  --wandb_project PointBAP --wandb_run crystal_regression_crossdocked_least_casf

#python point_vs/inference.py ../docked_regression/PointBAP/docked_regression/checkpoints/ckpt_epoch_1.pt /data/localhost/not-backed-up/durant/casf_2016_pymol_redocked.types /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets  --wandb_project PointBAP --wandb_run docked_regression_redocked_casf
#python point_vs/inference.py ../docked_regression/PointBAP/docked_regression/checkpoints/ckpt_epoch_1.pt /data/localhost/not-backed-up/durant/casf_2016_pymol_crossdocked_most_similiar_structure.types /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets  --wandb_project PointBAP --wandb_run docked_regression_crossdocked_most_casf
#python point_vs/inference.py ../docked_regression/PointBAP/docked_regression/checkpoints/ckpt_epoch_1.pt /data/localhost/not-backed-up/durant/casf_2016_pymol_crossdocked_least_similiar_structure.types /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets  --wandb_project PointBAP --wandb_run docked_regression_crossdocked_least_casf

python point_vs/inference_mass_run.py








echo The script has finised running
