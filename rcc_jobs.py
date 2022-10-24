import os
import math
phis = [0.01,0.05,0.1,0.2,0.5]
taus = [0.13]
maxdims = [200]
hz = 0.1
N_taus = [20]
numkicks = [400]
L = 100

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
sites_and_hzs = "["
for i in range(1,L+1):
    sites_and_hzs=str(sites_and_hzs)+"["+str(i)+","+str(hz)+"],"

sites_and_hzs = sites_and_hzs.rstrip(sites_and_hzs[-1])
sites_and_hzs=str(sites_and_hzs)+"]"

print("sites_and_hzs = ",sites_and_hzs)


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

                    
                    filename = "in_TEBD_open_staggered_delta_central_kicked_L_"+str(L)+"_hz_"+str(hz)+"_phi_"+str(phi)+"_tau_"+str(tau)+"_maxdim_"+str(maxdim)+"_numkicks_"+str(numkick)+"_Ntau_"+str(N_tau)
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
                    f.write("file_name=\"OPEN_TEBD_delta_central_kicked_L_"+str(L)+"_hz_"+str(hz)+"_phi_"+str(phi)+"_tau_"+str(tau)+"_maxdim_"+str(maxdim)+"_numkicks_"+str(numkick)+"_Ntau_"+str(N_tau)+"\"\n")
                    f.write("N_tau="+str(N_tau)+"\n")
                    f.close()
                    jobfname = "job"+str(filename) 
                    j = open(jobfname, 'w')
                    j.write("#!/bin/bash -l\n")
                    j.write("#Batch Queue Script\n")
                    j.write("#SBATCH --time=100:00:00\n")
                    j.write("#SBATCH --nodes=1 \n ")
                    j.write("#SBATCH --ntasks-per-node=1\n")
                    j.write("#SBATCH --cpus-per-task=16\n")
                    #j.write("#SBATCH -p genacc_q\n")
                    j.write("#SBATCH -p changlani_q\n")
                    j.write("export MKL_NUM_THREADS=16\n")
                    j.write("julia /gpfs/home/hchanglani/TEBDdrive/HJC_DELTA.jl "+str(filename)+" > "+"out_TEBD_open_staggered_delta_central_kicked_L_"+str(L)+"_hz_"+str(hz)+"_phi_"+str(phi)+"_tau_"+str(tau)+"_maxdim_"+str(maxdim)+"_numkicks_"+str(numkick)+"_Ntau_"+str(N_tau))
                    j.close()
                    os.system("chmod 755 "+jobfname)
                    os.system("sbatch "+jobfname)


