import os
import math

phis = [15]
taus = [0.834420]
maxdims = [100,200,400]
numkicks = [1000]
hzs = [1]
kicks = [1,5,9,13,17,21,25,29]


N_taus = [100]



dt = taus[0]/N_taus[0]
def checkdt(dt,Ntaus,tau):
    if dt > 0.01:
        N_taus = Ntau*10
        dt = tau/Ntau
        return checkdt(dt,Ntau)
    else:
        return dt,Ntaus

#dt,N_tau = checkdt(dt,N_tau,tau) 	



L = 30
sites_and_initial_product_states= "["
for i in range(1,L+1):
    sites_and_initial_product_states = str(sites_and_initial_product_states)+"["+str(i)+",[1/sqrt(2.0),0.0],[-1/sqrt(2.0),0.0]],"


sites_and_initial_product_states = sites_and_initial_product_states.rstrip(sites_and_initial_product_states[-1])
sites_and_initial_product_states = str(sites_and_initial_product_states)+"]"




print(" ")
print("sites_and_initial_product_states =",sites_and_initial_product_states)
print(" ")

bonds_and_Js= "["
for i in range(1,L):
    bonds_and_Js=str(bonds_and_Js)+"["+str(i)+","+str(i+1)+","+str((-1.0)**i)+","+str((-1.0)**i)+"],"
# set J = 0 for boundary for open
bonds_and_Js=str(bonds_and_Js)+"["+str(L)+","+str(1)+",0.0,0.0]]"

print("bonds_and_Js=",bonds_and_Js)
print(" ")
hz = 1
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
                        dt = tau/N_tau
                        dt,N_tau =  checkdt(dt,N_tau,tau)
                        for i in range(1,L+1):
                            k = 1
                            #sites_and_phis=str(sites_and_phis)+"["+str(ip)+",0.0]," 
                            for j in range(0,len(kicks)):
                                if i==kicks[j]:
                                    k = 0
                                    break
                            
                            if k == 1:
                                sites_and_phis=str(sites_and_phis)+"["+str(i)+",0.0],"
                                #break
                            if k == 0:
                                sites_and_phis=str(sites_and_phis)+"["+str(i)+","+str(phi)+"],"
                                #break

                        sites_and_phis = sites_and_phis.rstrip(sites_and_phis[-1])
                        sites_and_phis=str(sites_and_phis)+"]"
                        print(" ")
                        print("sites_and_phis = ",sites_and_phis)
#                        STOP
                        
                        filename = "file_Bhaskarian_pulse_pulse_40_site_"+str(kicks[0])+"_site_"+str(kicks[1])+"_site_"+str(kicks[2])+"_site_"+str(kicks[3])+"_site_"+str(kicks[4])+"_site_"+str(kicks[5])+"_site_"+str(kicks[6])+"_site_"+str(kicks[7])+"_L_"+str(L)+"_40_pulses_In_Bhaskarian_pulse_site_1_hz_"+str(hz)+"_phi_"+str(phi)+"_tau_"+str(tau)+"_maxdim_"+str(maxdim)+"_numkicks_"+str(numkick)+"_Ntau_"+str(N_tau)
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
                        f.write("file_name=\"OPEN_Bhaskarian_pulse_pulse_40_site_"+str(kicks[0])+"_site_"+str(kicks[1])+"_site_"+str(kicks[2])+"_site_"+str(kicks[3])+"_site_"+str(kicks[4])+"_site_"+str(kicks[5])+"_site_"+str(kicks[6])+"_site_"+str(kicks[7])+"_output_with_effective_ham_double_S_experiment_1_site_1_kicked_L_"+str(L)+"__hz_"+str(hz)+"_phi_"+str(phi)+"_tau_"+str(tau)+"_maxdim_"+str(maxdim)+"_numkicks_"+str(numkick)+"_Ntau_"+str(N_tau)+"\"\n")
                        f.write("N_tau="+str(N_tau)+"\n")
                        jobfname = "job"+str(filename) 
                        j = open(jobfname, 'w')
                        j.write("#!/bin/bash -l\n")
                        j.write("#Batch Queue Script\n")
                        j.write("#SBATCH --time=14-00:00:00\n")
                        j.write("#SBATCH --nodes=1 \n ")
                        j.write("#SBATCH --ntasks-per-node=1\n")
                        j.write("#SBATCH --cpus-per-task=16\n")
                        j.write("#SBATCH --mem=64G\n")
                        j.write("export OMP_NUM_THREADS=16 \n")
                        #j.write("export MKL_NUM_THREADS=16 \n")
                        j.write("#SBATCH -p genacc_q\n")
                        #j.write("#SBATCH -p changlani_q\n")
                        j.write("julia square_beta_tebd.jl "+str(filename)+" > "+"OPEN_pipe_out_OPEN_Bhaskarian_pulse_pulse_40_site_"+str(kicks[0])+"_site_"+str(kicks[1])+"_site_"+str(kicks[2])+"_site_"+str(kicks[3])+"_site_"+str(kicks[4])+"_site_"+str(kicks[5])+"_site_"+str(kicks[6])+"_site_"+str(kicks[7])+"_kicked_L_"+str(L)+"__hz_"+str(hz)+"_phi_"+str(phi)+"_tau_"+str(tau)+"_maxdim_"+str(maxdim)+"_numkicks_"+str(numkick)+"_Ntau_"+str(N_tau))
                        j.close()
                        os.system("chmod 755 "+jobfname)
                        os.system("sbatch "+jobfname)


