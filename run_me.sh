#!/bin/bash 
# MOAB/Torque submission script for SciNet GPC
#
#PBS -l nodes=1:ppn=8,walltime=02:00:00
#PBS -N lshape_ipcs_ab_cn_lshape_2_ts1000_cycles3_uOrder1

#PBS -q debug
# DIRECTORY TO RUN - $PBS_O_WORKDIR is directory job was submitted from
cd $PBS_O_WORKDIR
#export INSTANT_CACHE_DIR=Work/instant-cache/
#export INSTANT_ERROR_DIR=Work/instant-error/
# EXECUTION COMMAND; -np = nodes*ppn
# mpirun --tag-output -np 8 env  | grep LD_LIBRARY_PATH
mpirun -n 8 python ns lshape ipcs_ab_cn wct_hrs=2 tOrder=None period=951 node_type=debug restart_path=./results/lshape_ipcs_ab_cn_lshape_2_ts1000_cycles3_uOrder1 bc_out_type=zero_pressure bc_in_profile=parabolic bc_in_data=['FC_MCA'] processors=8 save_frequency=5 timesteps=1000 check_point_frequency=10 viscosity=0.0035 solver_name=ipcs_ab_cn case_name=lshape_ipcs_ab_cn_lshape_2_ts1000_cycles3_uOrder1 restart_time=101 wct_mins=0 restart=False problem_name=lshape mesh_name=lshape_2 no_of_restart=0 time_out=2850 cycles=3 uOrder=1  &> logs/lshape_ipcs_ab_cn_lshape_2_ts1000_cycles3_uOrder1_restart_0