#!/bin/bash   
#SBATCH -J regression                     # Job name
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
#python point_vs.py egnn /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets crystal_regression -b 16 -lr 0.0008 -e 10 --activation relu -k 32 --radius 6 --compact --weight_decay 1e-4 --train_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_train.types --min_inactive_rmsd 2 --max_active_rmsd 2 --egnn_residual --egnn_normalise --egnn_tanh --edge_radius 10 --val_on_epoch_end --norm_coords --norm_feats --warm_restarts --egnn_attention --dropout 0.0 --egnn_classify_on_feats --estimate_bonds --only_save_best_models --layers 48 --model_task multi_regression --load_weights /data/localhost/not-backed-up/durant/jack_models/sweep/48L_compact__0/checkpoints/ckpt_epoch_7.pt --test_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_val.types -t /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets --wandb_project PointBAP --wandb_run crystal_regression --top1

#python point_vs.py egnn /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets crystal_regression -b 16 -lr 0.0008 -e 10 --activation relu -k 32 --radius 6 --compact --weight_decay 1e-4 --train_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_sequence_threshold_100_train.types --min_inactive_rmsd 2 --max_active_rmsd 2 --egnn_residual --egnn_normalise --egnn_tanh --edge_radius 10 --val_on_epoch_end --norm_coords --norm_feats --warm_restarts --egnn_attention --dropout 0.0 --egnn_classify_on_feats --estimate_bonds --only_save_best_models --layers 48 --model_task multi_regression --load_weights /data/localhost/not-backed-up/durant/jack_models/sweep/48L_compact__0/checkpoints/ckpt_epoch_7.pt --test_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_sequence_threshold_100_val.types -t /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets --wandb_project PointBAP --wandb_run crystal_regression_seq_100 --top1

#python point_vs.py egnn /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets crystal_regression -b 16 -lr 0.0008 -e 10 --activation relu -k 32 --radius 6 --compact --weight_decay 1e-4 --train_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_sequence_threshold_90_train.types --min_inactive_rmsd 2 --max_active_rmsd 2 --egnn_residual --egnn_normalise --egnn_tanh --edge_radius 10 --val_on_epoch_end --norm_coords --norm_feats --warm_restarts --egnn_attention --dropout 0.0 --egnn_classify_on_feats --estimate_bonds --only_save_best_models --layers 48 --model_task multi_regression --load_weights /data/localhost/not-backed-up/durant/jack_models/sweep/48L_compact__0/checkpoints/ckpt_epoch_7.pt --test_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_sequence_threshold_90_val.types -t /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets --wandb_project PointBAP --wandb_run crystal_regression_seq_90 --top1

#python point_vs.py egnn /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets crystal_regression -b 16 -lr 0.0008 -e 10 --activation relu -k 32 --radius 6 --compact --weight_decay 1e-4 --train_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_sequence_threshold_30_train.types --min_inactive_rmsd 2 --max_active_rmsd 2 --egnn_residual --egnn_normalise --egnn_tanh --edge_radius 10 --val_on_epoch_end --norm_coords --norm_feats --warm_restarts --egnn_attention --dropout 0.0 --egnn_classify_on_feats --estimate_bonds --only_save_best_models --layers 48 --model_task multi_regression --load_weights /data/localhost/not-backed-up/durant/jack_models/sweep/48L_compact__0/checkpoints/ckpt_epoch_7.pt --test_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_sequence_threshold_30_val.types -t /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets --wandb_project PointBAP --wandb_run crystal_regression_seq_30 --top1

python point_vs.py egnn /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets crystal_regression -b 16 -lr 0.0008 -e 10 --activation relu -k 32 --radius 6 --compact --weight_decay 1e-4 --train_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_tanimoto_threshold_100_train.types --min_inactive_rmsd 2 --max_active_rmsd 2 --egnn_residual --egnn_normalise --egnn_tanh --edge_radius 10 --val_on_epoch_end --norm_coords --norm_feats --warm_restarts --egnn_attention --dropout 0.0 --egnn_classify_on_feats --estimate_bonds --only_save_best_models --layers 48 --model_task multi_regression --load_weights /data/localhost/not-backed-up/durant/jack_models/sweep/48L_compact__0/checkpoints/ckpt_epoch_7.pt --test_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_tanimoto_threshold_100_val.types -t /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets --wandb_project PointBAP --wandb_run crystal_regression_tan_100 --top1

#python point_vs.py egnn /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets crystal_regression -b 16 -lr 0.0008 -e 10 --activation relu -k 32 --radius 6 --compact --weight_decay 1e-4 --train_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_tanimoto_threshold_90_train.types --min_inactive_rmsd 2 --max_active_rmsd 2 --egnn_residual --egnn_normalise --egnn_tanh --edge_radius 10 --val_on_epoch_end --norm_coords --norm_feats --warm_restarts --egnn_attention --dropout 0.0 --egnn_classify_on_feats --estimate_bonds --only_save_best_models --layers 48 --model_task multi_regression --load_weights /data/localhost/not-backed-up/durant/jack_models/sweep/48L_compact__0/checkpoints/ckpt_epoch_7.pt --test_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_tanimoto_threshold_90_val.types -t /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets --wandb_project PointBAP --wandb_run crystal_regression_tan_90 --top1

#python point_vs.py egnn /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets crystal_regression -b 16 -lr 0.0008 -e 10 --activation relu -k 32 --radius 6 --compact --weight_decay 1e-4 --train_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_tanimoto_threshold_30_train.types --min_inactive_rmsd 2 --max_active_rmsd 2 --egnn_residual --egnn_normalise --egnn_tanh --edge_radius 10 --val_on_epoch_end --norm_coords --norm_feats --warm_restarts --egnn_attention --dropout 0.0 --egnn_classify_on_feats --estimate_bonds --only_save_best_models --layers 48 --model_task multi_regression --load_weights /data/localhost/not-backed-up/durant/jack_models/sweep/48L_compact__0/checkpoints/ckpt_epoch_7.pt --test_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_tanimoto_threshold_30_val.types -t /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets --wandb_project PointBAP --wandb_run crystal_regression_tan_30 --top1

#python point_vs.py egnn /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets crystal_regression -b 16 -lr 0.0008 -e 10 --activation relu -k 32 --radius 6 --compact --weight_decay 1e-4 --train_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_both_threshold_100_train.types --min_inactive_rmsd 2 --max_active_rmsd 2 --egnn_residual --egnn_normalise --egnn_tanh --edge_radius 10 --val_on_epoch_end --norm_coords --norm_feats --warm_restarts --egnn_attention --dropout 0.0 --egnn_classify_on_feats --estimate_bonds --only_save_best_models --layers 48 --model_task multi_regression --load_weights /data/localhost/not-backed-up/durant/jack_models/sweep/48L_compact__0/checkpoints/ckpt_epoch_7.pt --test_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_both_threshold_100_val.types -t /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets --wandb_project PointBAP --wandb_run crystal_regression_both_100 --top1

python point_vs.py egnn /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets crystal_regression -b 16 -lr 0.0008 -e 10 --activation relu -k 32 --radius 6 --compact --weight_decay 1e-4 --train_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_both_threshold_90_train.types --min_inactive_rmsd 2 --max_active_rmsd 2 --egnn_residual --egnn_normalise --egnn_tanh --edge_radius 10 --val_on_epoch_end --norm_coords --norm_feats --warm_restarts --egnn_attention --dropout 0.0 --egnn_classify_on_feats --estimate_bonds --only_save_best_models --layers 48 --model_task multi_regression --load_weights /data/localhost/not-backed-up/durant/jack_models/sweep/48L_compact__0/checkpoints/ckpt_epoch_7.pt --test_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_both_threshold_90_val.types -t /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets --wandb_project PointBAP --wandb_run crystal_regression_both_90 --top1

#python point_vs.py egnn /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets crystal_regression -b 16 -lr 0.0008 -e 10 --activation relu -k 32 --radius 6 --compact --weight_decay 1e-4 --train_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_both_threshold_30_train.types --min_inactive_rmsd 2 --max_active_rmsd 2 --egnn_residual --egnn_normalise --egnn_tanh --edge_radius 10 --val_on_epoch_end --norm_coords --norm_feats --warm_restarts --egnn_attention --dropout 0.0 --egnn_classify_on_feats --estimate_bonds --only_save_best_models --layers 48 --model_task multi_regression --load_weights /data/localhost/not-backed-up/durant/jack_models/sweep/48L_compact__0/checkpoints/ckpt_epoch_7.pt --test_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_both_threshold_30_val.types -t /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets --wandb_project PointBAP --wandb_run crystal_regression_both_30 --top1


#python point_vs.py egnn /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets crystal_regression -b 16 -lr 0.0008 -e 10 --activation relu -k 32 --radius 6 --compact --weight_decay 1e-4 --train_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_both_threshold_30_size_control_train.types --min_inactive_rmsd 2 --max_active_rmsd 2 --egnn_residual --egnn_normalise --egnn_tanh --edge_radius 10 --val_on_epoch_end --norm_coords --norm_feats --warm_restarts --egnn_attention --dropout 0.0 --egnn_classify_on_feats --estimate_bonds --only_save_best_models --layers 48 --model_task multi_regression --load_weights /data/localhost/not-backed-up/durant/jack_models/sweep/48L_compact__0/checkpoints/ckpt_epoch_7.pt --test_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_both_threshold_30_size_control_val.types -t /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets --wandb_project PointBAP --wandb_run crystal_regression_both_30_size_control --top1











echo Now running docking models


#python point_vs.py egnn /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets docked_regression -b 16 -lr 0.0008 -e 10 --activation relu -k 32 --radius 6 --compact --weight_decay 1e-4 --train_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_redocked_train.types --min_inactive_rmsd 2 --max_active_rmsd 2 --egnn_residual --egnn_normalise --egnn_tanh --edge_radius 10 --val_on_epoch_end --norm_coords --norm_feats --warm_restarts --egnn_attention --dropout 0.0 --egnn_classify_on_feats --estimate_bonds --only_save_best_models --layers 48 --model_task multi_regression --load_weights /data/localhost/not-backed-up/durant/jack_models/sweep/48L_compact__0/checkpoints/ckpt_epoch_7.pt --test_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_redocked_val.types -t /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets --wandb_project PointBAP --wandb_run docked_regression --top1

#python point_vs.py egnn /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets docked_regression -b 16 -lr 0.0008 -e 10 --activation relu -k 32 --radius 6 --compact --weight_decay 1e-4 --train_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_redocked_sequence_threshold_100_train.types --min_inactive_rmsd 2 --max_active_rmsd 2 --egnn_residual --egnn_normalise --egnn_tanh --edge_radius 10 --val_on_epoch_end --norm_coords --norm_feats --warm_restarts --egnn_attention --dropout 0.0 --egnn_classify_on_feats --estimate_bonds --only_save_best_models --layers 48 --model_task multi_regression --load_weights /data/localhost/not-backed-up/durant/jack_models/sweep/48L_compact__0/checkpoints/ckpt_epoch_7.pt --test_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_redocked_sequence_threshold_100_val.types  -t /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets --wandb_project PointBAP --wandb_run docked_regression_seq_100 --top1

#python point_vs.py egnn /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets docked_regression -b 16 -lr 0.0008 -e 10 --activation relu -k 32 --radius 6 --compact --weight_decay 1e-4 --train_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_redocked_sequence_threshold_90_train.types --min_inactive_rmsd 2 --max_active_rmsd 2 --egnn_residual --egnn_normalise --egnn_tanh --edge_radius 10 --val_on_epoch_end --norm_coords --norm_feats --warm_restarts --egnn_attention --dropout 0.0 --egnn_classify_on_feats --estimate_bonds --only_save_best_models --layers 48 --model_task multi_regression --load_weights /data/localhost/not-backed-up/durant/jack_models/sweep/48L_compact__0/checkpoints/ckpt_epoch_7.pt --test_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_redocked_sequence_threshold_90_val.types -t /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets --wandb_project PointBAP --wandb_run docked_regression_seq_90 --top1

python point_vs.py egnn /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets docked_regression -b 16 -lr 0.0008 -e 10 --activation relu -k 32 --radius 6 --compact --weight_decay 1e-4 --train_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_redocked_sequence_threshold_30_train.types --min_inactive_rmsd 2 --max_active_rmsd 2 --egnn_residual --egnn_normalise --egnn_tanh --edge_radius 10 --val_on_epoch_end --norm_coords --norm_feats --warm_restarts --egnn_attention --dropout 0.0 --egnn_classify_on_feats --estimate_bonds --only_save_best_models --layers 48 --model_task multi_regression --load_weights /data/localhost/not-backed-up/durant/jack_models/sweep/48L_compact__0/checkpoints/ckpt_epoch_7.pt --test_types  /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_redocked_sequence_threshold_30_val.types -t /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets --wandb_project PointBAP --wandb_run docked_regression_seq_30 --top1

#python point_vs.py egnn /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets docked_regression -b 16 -lr 0.0008 -e 10 --activation relu -k 32 --radius 6 --compact --weight_decay 1e-4 --train_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_redocked_tanimoto_threshold_100_train.types --min_inactive_rmsd 2 --max_active_rmsd 2 --egnn_residual --egnn_normalise --egnn_tanh --edge_radius 10 --val_on_epoch_end --norm_coords --norm_feats --warm_restarts --egnn_attention --dropout 0.0 --egnn_classify_on_feats --estimate_bonds --only_save_best_models --layers 48 --model_task multi_regression --load_weights /data/localhost/not-backed-up/durant/jack_models/sweep/48L_compact__0/checkpoints/ckpt_epoch_7.pt --test_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_redocked_tanimoto_threshold_100_val.types -t /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets --wandb_project PointBAP --wandb_run docked_regression_tan_100 --top1

#python point_vs.py egnn /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets docked_regression -b 16 -lr 0.0008 -e 10 --activation relu -k 32 --radius 6 --compact --weight_decay 1e-4 --train_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_redocked_tanimoto_threshold_90_train.types --min_inactive_rmsd 2 --max_active_rmsd 2 --egnn_residual --egnn_normalise --egnn_tanh --edge_radius 10 --val_on_epoch_end --norm_coords --norm_feats --warm_restarts --egnn_attention --dropout 0.0 --egnn_classify_on_feats --estimate_bonds --only_save_best_models --layers 48 --model_task multi_regression --load_weights /data/localhost/not-backed-up/durant/jack_models/sweep/48L_compact__0/checkpoints/ckpt_epoch_7.pt --test_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_redocked_tanimoto_threshold_90_val.types -t /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets --wandb_project PointBAP --wandb_run docked_regression_tan_90 --top1

#python point_vs.py egnn /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets docked_regression -b 16 -lr 0.0008 -e 10 --activation relu -k 32 --radius 6 --compact --weight_decay 1e-4 --train_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_redocked_tanimoto_threshold_30_train.types --min_inactive_rmsd 2 --max_active_rmsd 2 --egnn_residual --egnn_normalise --egnn_tanh --edge_radius 10 --val_on_epoch_end --norm_coords --norm_feats --warm_restarts --egnn_attention --dropout 0.0 --egnn_classify_on_feats --estimate_bonds --only_save_best_models --layers 48 --model_task multi_regression --load_weights /data/localhost/not-backed-up/durant/jack_models/sweep/48L_compact__0/checkpoints/ckpt_epoch_7.pt --test_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_redocked_tanimoto_threshold_30_val.types -t /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets --wandb_project PointBAP --wandb_run docked_regression_tan_30 --top1

#python point_vs.py egnn /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets docked_regression -b 16 -lr 0.0008 -e 10 --activation relu -k 32 --radius 6 --compact --weight_decay 1e-4 --train_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_redocked_both_threshold_100_train.types --min_inactive_rmsd 2 --max_active_rmsd 2 --egnn_residual --egnn_normalise --egnn_tanh --edge_radius 10 --val_on_epoch_end --norm_coords --norm_feats --warm_restarts --egnn_attention --dropout 0.0 --egnn_classify_on_feats --estimate_bonds --only_save_best_models --layers 48 --model_task multi_regression --load_weights /data/localhost/not-backed-up/durant/jack_models/sweep/48L_compact__0/checkpoints/ckpt_epoch_7.pt --test_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_redocked_both_threshold_100_val.types -t /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets --wandb_project PointBAP --wandb_run docked_regression_both_100 --top1

#python point_vs.py egnn /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets docked_regression -b 16 -lr 0.0008 -e 10 --activation relu -k 32 --radius 6 --compact --weight_decay 1e-4 --train_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_redocked_both_threshold_90_train.types --min_inactive_rmsd 2 --max_active_rmsd 2 --egnn_residual --egnn_normalise --egnn_tanh --edge_radius 10 --val_on_epoch_end --norm_coords --norm_feats --warm_restarts --egnn_attention --dropout 0.0 --egnn_classify_on_feats --estimate_bonds --only_save_best_models --layers 48 --model_task multi_regression --load_weights /data/localhost/not-backed-up/durant/jack_models/sweep/48L_compact__0/checkpoints/ckpt_epoch_7.pt --test_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_redocked_both_threshold_90_val.types -t /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets --wandb_project PointBAP --wandb_run docked_regression_both_90 --top1

#python point_vs.py egnn /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets docked_regression -b 16 -lr 0.0008 -e 10 --activation relu -k 32 --radius 6 --compact --weight_decay 1e-4 --train_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_redocked_both_threshold_30_train.types --min_inactive_rmsd 2 --max_active_rmsd 2 --egnn_residual --egnn_normalise --egnn_tanh --edge_radius 10 --val_on_epoch_end --norm_coords --norm_feats --warm_restarts --egnn_attention --dropout 0.0 --egnn_classify_on_feats --estimate_bonds --only_save_best_models --layers 48 --model_task multi_regression --load_weights /data/localhost/not-backed-up/durant/jack_models/sweep/48L_compact__0/checkpoints/ckpt_epoch_7.pt --test_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_redocked_both_threshold_30_val.types -t /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets --wandb_project PointBAP --wandb_run docked_regression_both_30 --top1

#python point_vs.py egnn /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets docked_regression -b 16 -lr 0.0008 -e 10 --activation relu -k 32 --radius 6 --compact --weight_decay 1e-4 --train_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_redocked_both_threshold_30_size_control_train.types --min_inactive_rmsd 2 --max_active_rmsd 2 --egnn_residual --egnn_normalise --egnn_tanh --edge_radius 10 --val_on_epoch_end --norm_coords --norm_feats --warm_restarts --egnn_attention --dropout 0.0 --egnn_classify_on_feats --estimate_bonds --only_save_best_models --layers 48 --model_task multi_regression --load_weights /data/localhost/not-backed-up/durant/jack_models/sweep/48L_compact__0/checkpoints/ckpt_epoch_7.pt --test_types /data/localhost/not-backed-up/durant/pdbbind_2020_general_pymol_redocked_both_threshold_30_size_control_val.types -t /data/localhost/not-backed-up/durant/pdbbind_2020_general_parquets --wandb_project PointBAP --wandb_run docked_regression_both_30_size_control --top1
echo The script has finised running