#! /bin/bash
# for data

# source /global/home/users/yuanjing/.bashrc
# conda activate env_root
# input=$1


input="H3L_KF.list"
# root -b -q readtree.C\(\"$input\",0,\"fout_H3L_data_KF_0_80.root\",0,0,8\)
# root -b -q readtree.C\(\"$input\",0,\"fout_H3L_data_KF_0_40.root\",0,4,8\)
root -b -q readtree.C\(\"$input\",0,\"fout_H3L_data_KF_0_10.root\",0,7,8\)
#
# input="H3L_KF.list"
# root -b -q readtree.C\(\"$input\",0,\"fout_H3L_data_KF_10_40.root\",0,4,6\)

input="H3L2b_rt.root"
# root -b -q readtree.C\(\"$input\",0,\"fout_H3L_data_RT_0_80.root\",0,0,8\)
# root -b -q readtree.C\(\"$input\",0,\"fout_H3L_data_RT_0_40.root\",0,4,8\)
root -b -q readtree.C\(\"$input\",0,\"fout_H3L_data_RT_0_10.root\",0,7,8\)
#
# input="H3L2b_rt.root"
# root -b -q readtree.C\(\"$input\",0,\"fout_H3L_data_RT_10_40.root\",0,4,6\)
#
#
input="H3L2b_mc.root"
# input="H3L2b_mc_262.root"
# root -b -q readMc.C\(\"$input\",0,\"fout_H3L_data_MC_0080.root\",1,0,8\)
# root -b -q readtree.C\(\"$input\",0,\"fout_H3L_data_MC_RC_0080.root\",1,0,8\)
# #
# root -b -q readMc.C\(\"$input\",0,\"fout_H3L_data_MC_0_40.root\",1,4,8\)
# root  readMc.C\(\"$input\",0,\"fout_H3L_data_MC_0_40.root\",1,4,8\)
# root -b -q readtree.C\(\"$input\",0,\"fout_H3L_data_MC_RC_0_40.root\",1,4,8\)
root -b -q readMc.C\(\"$input\",0,\"fout_H3L_data_MC_0_10.root\",1,7,8\)
root -b -q readtree.C\(\"$input\",0,\"fout_H3L_data_MC_RC_0_10.root\",1,7,8\)
#
#

