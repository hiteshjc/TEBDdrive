import os
import math
#taus = [11,12,13,14,15,16,17,18,19,21,22,23,24,25,26,27,28,29]

#phis = [0.00001,0.0001,0.001,0.01,0.1,1,10]
#phis = [0.01]
#phis = [10.000006]

#taus = [0.001]
#taus = [30]
#taus = [0.011,0.012,0.013,0.014,0.015,0.016,0.017,0.018,0.019] 
#taus = [0.8]

phis = [0.5]
taus = [0.1]
maxdims = [400]

#taus = [0.001,0.2]

#taus = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09]
#taus = [0.1,0.01,0.001,1,2,3,4,5,6,7,8,9]

#phis = [0.5]
#taus = [0.001]

#taus = [0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
#phis=[0.0001,0.001,0.01]

#taus = [0.000001]
#phis = [0.001]
#phis = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
#phis = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
#phis = [0.00001]
#phis = [0.0001,0.001,0.01,0.1] 
#phis = [0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09]

hzs = [0.1]
N_taus = [10]
numkicks = [400]


L = 60
sites_and_initial_product_states= "["
for i in range(1,L+1):
    sites_and_initial_product_states = str(sites_and_initial_product_states)+"["+str(i)+",[1/sqrt(2.0),0.0],[1/sqrt(2.0),0.0]],"


sites_and_initial_product_states = sites_and_initial_product_states.rstrip(sites_and_initial_product_states[-1])
sites_and_initial_product_states = str(sites_and_initial_product_states)+"]"




print(" ")
print("sites_and_initial_product_states =",sites_and_initial_product_states)
print(" ")

bonds_and_Js= "["
for i in range(1,L):
    bonds_and_Js=str(bonds_and_Js)+"["+str(i)+","+str(i+1)+","+str((-1.0)**i)+","+str((-1.0)**i)+"],"
bonds_and_Js=str(bonds_and_Js)+"["+str(L)+","+str(1)+","+str((-1.0)**L)+","+str((-1.0)**L)+"]]"

print("bonds_and_Js=",bonds_and_Js)
print(" ")
hz = 0.1
sites_and_hzs = "["
for i in range(1,L+1):
    sites_and_hzs=str(sites_and_hzs)+"["+str(i)+","+str(hz)+"],"

sites_and_hzs = sites_and_hzs.rstrip(sites_and_hzs[-1])


sites_and_hzs=str(sites_and_hzs)+"]"

print("sites_and_hzs = ",sites_and_hzs)


for hz in hzs:
    for tau in taus:
        for phi in phis:
            for maxdim in maxdims:
                for N_tau in N_taus:
                    for numkick in numkicks:
                        sites_and_phis ="["

                        for i in range(1,L+1):
                            sites_and_phis=str(sites_and_phis)+"["+str(i)+",0.0],"
                            # half chain site kick
                            if i == int(L/2):
                                sites_and_phis=str(sites_and_phis)+"["+str(i)+","+str(phi)+"],"
     
                        sites_and_phis = sites_and_phis.rstrip(sites_and_phis[-1])
                        sites_and_phis=str(sites_and_phis)+"]"
                        print(" ")
                        print("sites_and_phis = ",sites_and_phis)

                        
                        filename = "L_"+str(L)+"_In_OPEN_TEBD_S_i_S_j_staggered_slight_translation_breaking_with_effective_ham_double_S_experiment_site_1_kicked_hz_"+str(hz)+"_phi_"+str(phi)+"_tau_"+str(tau)+"_maxdim_"+str(maxdim)+"_numkicks_"+str(numkick)+"_Ntau_"+str(N_tau)
                        f = open(filename, 'w')
                        f.write("L="+str(L)+"\n")
                        f.write("sites_and_initial_product_states="+str(sites_and_initial_product_states)+"\n")
                        f.write("bonds_and_Js="+str(bonds_and_Js)+"\n")
                        f.write("sites_and_hzs="+str(sites_and_hzs)+"\n")
                        f.write("tau="+str(tau)+"\n")
                        # f.write("sites_and_phis=[[1,"+str(phi)+"],[2,"+str(phi)+"],[3,"+str(phi)+"],[4,"+str(phi)+"],[5,"+str(phi)+"],[6,"+str(phi)+"],[7,"+str(phi)+"],[8,"+str(phi)+"],[9,"+str(phi)+"],[10,"+str(phi)+"]]\n")
                        f.write("sites_and_phis="+str(sites_and_phis)+"\n")
                        f.write("numkicks="+str(numkick)+"\n")
                        f.write("maxdim="+str(maxdim)+"\n")
                        f.write("cutoff=10^(-8)\n")
                        f.write("file_name=\"OPEN_TEBD_S_i_S_j_staggered_slight_translation_breaking_output_with_effective_ham_double_S_experiment_1_site_1_kicked_L_"+str(L)+"__hz_"+str(hz)+"_phi_"+str(phi)+"_tau_"+str(tau)+"_maxdim_"+str(maxdim)+"_numkicks_"+str(numkick)+"_Ntau_"+str(N_tau)+"\"\n")
                        f.write("N_tau="+str(N_tau)+"\n")
                        f.close()
                        #jobfname = "job"+str(filename) 
                        #j = open(jobfname, 'w')
                        #j.write("#!/bin/bash -l\n")
                        #j.write("#Batch Queue Script\n")
                        #j.write("#SBATCH --time=14-00:00:00\n")
                        #j.write("#SBATCH --nodes=1 \n ")
                        #j.write("#SBATCH --ntasks-per-node=1\n")
                        #j.write("#SBATCH --cpus-per-task=16\n")
                        #j.write("#SBATCH --mem=64G\n")
                      #  j.write("export OMP_NUM_THREADS=8 \n")
                      #  j.write("export MKL_NUM_THREADS=24 \n")
                        #j.write("#SBATCH -p genacc_q\n")
                        #j.write("#SBATCH -p changlani_q\n")
                        os.system("/home/rmelendrez/julia/julia TEBD_HPC_Staggered_Heisen.jl "+str(filename)+" > "+"pipe_out_TEBD_S_i_S_j_staggered_slight_translation_breaking_output_with_effective_ham_double_S_experiment_1_site_1_kicked_L_"+str(L)+"__hz_"+str(hz)+"_phi_"+str(phi)+"_tau_"+str(tau)+"_maxdim_"+str(maxdim)+"_numkicks_"+str(numkick)+"_Ntau_"+str(N_tau))
                        #os.system("chmod 755 "+jobfname)
                        #os.system("sbatch "+jobfname)


