# commented out XX corr
# Ronald Melendrez
# 1/23/2023
# SQUARE PULSE TEBD

include("single_site_entropy.jl")
include("Bhaskarian_gates_OBC.jl")

using ITensors
using JLD2
#using TimerOutputs


# function that reads in input data and stores it
function read_input_file(fname)
#
	data = []
	f = open(fname,"r")
#
	for ln in eachline(f)
		
		line = split(ln ,"=")
		push!(data,[line[1],Meta.parse(line[2])])

	end
	
	close(f)
	return data
end

# file name passed as an argument for this script
fname = ARGS[1]

# runs function and stores input data as a string literal inside of filedata
filedata = read_input_file(fname)






#############################################################
# Input Parameters
#############################################################

for value in filedata
	
	if (value[1]=="N_tau")	
	global Ntau = eval(value[2])
	end


	if(value[1]=="L") 
	global	N = eval(value[2])
	end
	
	if(value[1]=="sites_and_initial_product_states") 
	global	sites_and_initial_product_states = eval(value[2])
	end

	if(value[1]=="bonds_and_Js") 
	global	bonds_and_Js = eval(value[2])
	end
	
	if(value[1]=="maxdim") 
	global	maxdim1 = eval(value[2])
	end
	
	if(value[1]=="cutoff") 
	global	cutoff1 = eval(value[2])
	end
	
	if(value[1]=="sites_and_hzs") 
	global	sites_and_hzs = eval(value[2])
	end

	if(value[1]=="sites_and_phis") 
	global	sites_and_phis = eval(value[2])
	end
	
	if(value[1]=="tau") 
	global	tau = eval(value[2])
	end

	
	if(value[1]=="numkicks")
	global numkicks = eval(value[2])
	end
	
	if(value[1]=="file_name")
	global file_name = String(eval(value[2]))
	end

end





#data load complete statement
println("data load complete for ",file_name)


#################################################
#
# Start Extra Parameters
#
#################################################
#const to = TimerOutput()
s = siteinds("S=1/2",N;conserve_qns=false)




#calculating parameters from given parameters
println("Calculating parameters from input file...")
# using user inputs to calculate dt

# Rescaled period of pulse
tau = tau
#println("τ rescaled by 2π... ")

dt = tau/Ntau
println("calculated TEBD time evolution time step dt = ",dt)

# using user inputs to calculate ttotal
ttotal = numkicks*tau
println("Total time evolution T = ",ttotal)


#################################################
#
# End Extra Parameters
#
#################################################


# CALL GATE FUNCTION Located in Bhaskarian_gates.jl
gates1,gates2 = square_pulse_gate_generator(s,N,bonds_and_Js,sites_and_phis,sites_and_hzs,dt)


########################################################################
#
# Start State Construction
#
#######################################################################


states = []

    for n in 1:N
        push!(states, "Up")
    end

global psi0 = MPS(s,states)

println("Bond Dimension = ", linkdims(psi0))



###############################################
#
#  Polarize spins along initial directions
#
###############################################

for site_and_initial_product_state in sites_and_initial_product_states

        cosa =site_and_initial_product_state[2][1]
        sina =site_and_initial_product_state[3][1]

        α = site_and_initial_product_state[3][2]
        θ = angle(cosa+im*sina)
        i = site_and_initial_product_state[1]

        OP1 =cos(θ)*op("Id", s[i]) -2*im*sin(θ)*op("Sy",s[i])
        OP2 =(1+cis(α))/2 * op("Id", s[i]) +(1-cis(α))/2 * 2 * op("Sz", s[i])

        global psi0 = apply(OP1,psi0;cutoff1,maxdim1)
        normalize(psi0)
        global psi0 = apply(OP2,psi0;cutoff1,maxdim1)
        normalize(psi0)
end


# make two MPSs one which is all X polarized state
# another which will be time evolved and one which will
# remain |ψ(0)>
# here |ψ₁> ≡ exp(-i H₂τ/2) exp(-i H₁τ/2) |ψ(0)> = |ψ(τ)>
# and  |ψ₂> = |ψ(0)>

psi1 = deepcopy(psi0)
psi2 = deepcopy(psi1)

# initialize data arrays
# for time series


#######################################################
#
# End of State Construction
#
######################################################


#####################################################
#
# Start of Time Dynamics
#
#####################################################



println("initializing time series arrays!")
time = []
Prob = []
Bonddim =[]
XXmat = []
Xmat = []
entropies = []
entropys = []

for σ in 1:N
    # entropy_von_neumann function located in single_site_entropy.jl
    S1 = entropy_von_neumann(psi1,σ)
    println("S",σ,"=",S1)
    push!(entropys,S1)

end

push!(entropies, entropys)

# set t = 0
global t = 0.0

# record Echo for t = 0
Echo = abs2(inner(psi1,psi2))

# record Magnization
magz = expect(psi1,"Sx")

# record xx correlation matrix
xxcorr= correlation_matrix(psi1,"Sx","Sx")

push!(XXmat,xxcorr)
push!(Bonddim, maxlinkdim(psi1))
push!(Xmat,magz)
push!(time,t)
push!(Prob,Echo)

println("initializing time series arrays!")

# We chop up τ into τ/N pieces which we call Ntau
# We want to time evolve
# by gates1
# from 0 to τ/2
# and gates2
# from τ/2 to τ

# moreover we want to take measurements every τ
# suppose N = 10
# need to time evolve with gates1 for n_tau = [1,2,3,4,5]
# and need to evolve with gates2 for n_tau = [6,7,8,9,10]

# initialize n_tau = 1
n_tau = 1
switch = true
println("hello world")
println("Ntau=",Ntau)
Ntotal = Ntau*numkicks
println("Ntotal=", Ntotal)
# Rinse repeat and time evolve
#flush(STDOUT)
for n in 1:Ntotal


        println("n=",n)
	flush(stdout)

        global t =dt+ t
        t≈ttotal && break

		#evolve with gates1 for 0 >=  t <= τ/2
                if n_tau <= Ntau/2
                        global psi1 = apply(gates1, psi1; cutoff1,maxdim1)
                        psi1=normalize(psi1)

                end

		#evolve with gates2 for τ >= t > τ/2
                if n_tau > Ntau/2
                        global psi1 = apply(gates2, psi1; cutoff1,maxdim1)
                        psi1 = normalize(psi1)

                        if n_tau == Ntau
                                
					#Stroboscopic Time Append
					#push!(time,t)
                                        
					#Measure XX correlations
				#	xxcorr= correlation_matrix(psi1,"Sx","Sx")
                                 #       push!(XXmat,xxcorr)
                                        
					#Measure on site Magnitization
					magz = expect(psi1,"Sx")
                                    	#push!(Xmat,magz)
					
					#Measure fidelity
                                        local Echo = abs2(inner(psi1,psi2))
                                        #push!(Prob,Echo)
                                        
					#Measure Max bond dim accumulated
					#push!(Bonddim, maxlinkdim(psi1))
                                        println("t = ","maxbonddim = ", maxlinkdim(psi1))
					println(" ") 	

					
					#Measure on site entropy 
                                        local entropys = []
                                        
					for σ in 1:N
                                        #        push!(entropys , entropy_von_neumann(psi1,σ))
                                        end
					
					#push!(entropies, entropys)
                                        
					# reinitialize n_tau
					n_tau = 1

                            global switch = false
                        end

                end



                if switch ==  true
                        global  n_tau = 1 + n_tau
                end

                switch = true


    end

##########################################################################
#
# End of Time Dynamics
#
#########################################################################

time = time


#file_name = "data_1"

#for ele in Bonddim

#println("Bonddim = ", ele)

#end


@save file_name time Prob Bonddim Xmat XXmat entropies
