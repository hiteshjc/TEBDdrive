import os
import math
phis = [0.0000001]
taus = [0.001]
maxdims = [32]

hzs = [0.1]
N_taus = [10]
numkicks = [40000]


L = 10


#############################################
# Sets up initial product |X> state in loop
#############################################


sites_and_initial_product_states = "["
for i in range(1, L+1):
    sites_and_initial_product_states = str(
        sites_and_initial_product_states)+"["+str(i)+",[1/sqrt(2.0),0.0],[1/sqrt(2.0),0.0]],"


sites_and_initial_product_states = sites_and_initial_product_states.rstrip(
    sites_and_initial_product_states[-1])
sites_and_initial_product_states = str(sites_and_initial_product_states)+"]"
###############################################


#################################################################
# prints initial product states as will be produced in input file
#################################################################
print(" ")
print("sites_and_initial_product_states =", sites_and_initial_product_states)
print(" ")

##################################################################
# Sets up bonds "staggered" in loop
##################################################################
bonds_and_Js = "["
for i in range(1, L):
    bonds_and_Js = str(
        bonds_and_Js)+"["+str(i)+","+str(i+1)+","+str((-1.0)**i)+","+str((-1.0)**i)+"],"
bonds_and_Js = str(
    bonds_and_Js)+"["+str(L)+","+str(1)+","+str((-1.0)**L)+","+str((-1.0)**L)+"]]"

###############################################################
# Prints bonds as will be produced in input file
###############################################################

print("bonds_and_Js=", bonds_and_Js)
print(" ")


#################################################################
# print statement for how hz field will be produced in input file
#################################################################
hz = 0.1
sites_and_hzs = "["
for i in range(1, L+1):
    sites_and_hzs = str(sites_and_hzs)+"["+str(i)+","+str(hz)+"],"

sites_and_hzs = sites_and_hzs.rstrip(sites_and_hzs[-1])


sites_and_hzs = str(sites_and_hzs)+"]"

print("sites_and_hzs = ", sites_and_hzs)

#################################################################
# loop structure begins to generate input files here
#################################################################
for hz in hzs:
    for tau in taus:
        for phi in phis:
            for maxdim in maxdims:
                for N_tau in N_taus:
                    for numkick in numkicks:

                        ##################################################
                        #  Beginning drive protocol implemented input
                        ##################################################

                        sites_and_phis = "["

                        for i in range(1, L+1):
                            sites_and_phis = str(
                                sites_and_phis)+"["+str(i)+",0.0],"
                            # half chain site kick
                            if i == int(L/2):
                                sites_and_phis = str(
                                    sites_and_phis)+"["+str(i)+","+str(phi)+"],"

                        sites_and_phis = sites_and_phis.rstrip(
                            sites_and_phis[-1])
                        sites_and_phis = str(sites_and_phis)+"]"

                        ##################################################
                        # Printing drive protocol input
                        #################################################

                        print(" ")
                        print("sites_and_phis = ", sites_and_phis)

                        filename = "L_"+str(L)+"_In_TEBD_SQUARE_with_SVN_experiment_site_1_kicked_hz_"+str(hz)+"_phi_"+str(
                            phi)+"_tau_"+str(tau)+"_maxdim_"+str(maxdim)+"_numkicks_"+str(numkick)+"_Ntau_"+str(N_tau)
                        f = open(filename, 'w')
                        f.write("L="+str(L)+"\n")
                        f.write("sites_and_initial_product_states=" +
                                str(sites_and_initial_product_states)+"\n")
                        f.write("bonds_and_Js="+str(bonds_and_Js)+"\n")
                        f.write("sites_and_hzs="+str(sites_and_hzs)+"\n")
                        f.write("tau="+str(tau)+"\n")
                        # f.write("sites_and_phis=[[1,"+str(phi)+"],[2,"+str(phi)+"],[3,"+str(phi)+"],[4,"+str(phi)+"],[5,"+str(phi)+"],[6,"+str(phi)+"],[7,"+str(phi)+"],[8,"+str(phi)+"],[9,"+str(phi)+"],[10,"+str(phi)+"]]\n")
                        f.write("sites_and_phis="+str(sites_and_phis)+"\n")
                        f.write("numkicks="+str(numkick)+"\n")
                        f.write("maxdim="+str(maxdim)+"\n")
                        f.write("cutoff=10^(-8)\n")
                        f.write("file_name=\"TEBD_SQUARE_with_SVN_output_experiment_1_site_1_kicked_L_"+str(L)+"__hz_"+str(hz)+"_phi_"+str(
                            phi)+"_tau_"+str(tau)+"_maxdim_"+str(maxdim)+"_numkicks_"+str(numkick)+"_Ntau_"+str(N_tau)+"\"\n")
                        f.write("N_tau="+str(N_tau)+"\n")
                        f.close()

                        # should install lastest version ofjulia locally on planck with version of ITensors
                        os.system("/home/rmelendrez/julia/julia entropy_hack_SQUARE_TEBD_HPC_Staggered_Heisen.jl "+str(filename)+" > "+"pipe_out_SQUARE_SVN_TEBD_output_experiment_1_site_1_kicked_L_"+str(
                            L)+"_hz_"+str(hz)+"_phi_"+str(phi)+"_tau_"+str(tau)+"_maxdim_"+str(maxdim)+"_numkicks_"+str(numkick)+"_Ntau_"+str(N_tau))
