#!/bin/bash   
#SBATCH -J redock                     # Job name
#SBATCH -A durant                         # Project Account                                   
#S BATCH --time=00:30:00                 # Walltime                                      
#SBATCH --mem-per-cpu=4096             # memory/cpu (in MB) ### commented out              
#SBATCH --ntasks=1                      # 1 tasks                                               
#SBATCH --cpus-per-task=8               # number of cores per task                          
#SBATCH --nodes=1 # number of nodes
#S BATCH --gres=gpu:1
#SBATCH --overcommit
# S BATCH --gpus-per-node=1
#S BATCH --exclusive                     # node should not be shared with other jobs, only use this if you intend the node to be usable only by you as this will block other users from submitting jobs to the same node                
#SBATCH --chdir=/data/localhost/not-backed-up/durant  # From where you want the job to be run
#SBATCH --mail-user=guy.durant@linacre.ox.ac.uk  # set email address                           
#S BATCH --mail-type=ALL                 # Spam us with everything, caution
#SBATCH --mail-type=begin               # Instead only email when job begins...
#SBATCH --mail-type=end                 # ... and ends
#SBATCH --clusters=srf_cpu_01
#SBATCH --partition=high-opig-cpu    # Select a specific partition rather than default 
#SBATCH -w naga02.cpu.stats.ox.ac.uk # Provide a specific node/nodelist rather than the standard nodelist associated with the partition (useful if you have a data setup on one specific node)
#SBATCH --output=/data/localhost/not-backed-up/durant/slurm_%j.out  # Writes standard output to this file. %j is jobnumber                             
#SBATCH --error=/data/localhost/not-backed-up/durant/slurm_%j.err   # Writes error messages to this file. %j is jobnumber

source ~/.bashrc
source ~/miniconda3/etc/profile.d/conda.sh
conda activate vina
echo The batch script has been starteD
#python3 point_vs.py egnn /data/localhost/not-backed-up/durant/CrossDocked2020_parquets/ test_output --train_types /data/localhost/not-backed-up/durant/redocked_types_train.types --test_types /data/localhost/not-backed-up/durant/redocked_types_val.types --egnn_classify_on_feats --epochs 100  --wandb_project PointBAP --wandb_run test_class --batch_size 32 --compact --use_1cycle
#python point_vs.py egnn /data/localhost/not-backed-up/durant/CrossDocked2020_parquets/ test_output -t /data/localhost/not-backed-up/durant/CrossDocked2020_parquets/ -b 64 -lr 0.001 -e 20 --activation relu -k 32 --radius 6 --wandb_project PointBAP --compact --weight_decay 1e-4 --train_types /data/localhost/not-backed-up/durant/redocked_types_train.types --test_types /data/localhost/not-backed-up/durant/redocked_types_val.types  --egnn_residual --egnn_normalise --egnn_tanh --edge_radius 10 --norm_coords --norm_feats  --warm_restarts --egnn_attention --dropout 0.0 --egnn_classify_on_feats --estimate_bonds  --layers 12  --wandb_run test_class --val_on_epoch_end
#python point_vs.py egnn /data/localhost/not-backed-up/durant/CrossDocked2020_parquets/ test_output -t /data/localhost/not-backed-up/durant/CrossDocked2020_parquets/ -b 32 -lr 0.001 -e 20 --activation relu -k 32 --radius 6 --wandb_project PointBAP --compact --weight_decay 1e-4 --train_types /data/localhost/not-backed-up/durant/redocked_types_train.types --test_types /data/localhost/not-backed-up/durant/redocked_types_val.types  --egnn_residual --egnn_normalise --egnn_tanh --edge_radius 10 --norm_coords --norm_feats  --warm_restarts --egnn_attention --dropout 0.0 --egnn_classify_on_feats --estimate_bonds  --layers 48  --wandb_run test_class_more_layers --val_on_epoch_end
python3 run_crossdocking_casf_vina.py
echo The script has finised running
