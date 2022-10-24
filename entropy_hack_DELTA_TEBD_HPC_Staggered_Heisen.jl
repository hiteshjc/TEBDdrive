using ITensors
using JLD2
using TimerOutputs
#using PyPlot
#using LaTeXStrings
using PrettyTables

const to = TimerOutput()

include("single_site_entropy.jl")

####################################################
####################################################

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



# Attributes variables that we want to get from the input file 
# Note that if you need to add additional variables
# you MUST update the attributes and extraction
# portion of script

attributes = ["L","N_tau","sites_and_inital_product_states","bonds_and_Js","sites_and_hzs","sites_and_phis","tau","numkicks","file_name"]

# file name passed as an argument for this script
fname = ARGS[1]

# runs function and stores input data as a string literal inside of filedata
filedata = read_input_file(fname)

	
#############################################################################
#############################################################################
#
# EXTRACTION OF SCRIPT DATA AND ASSIGNED TO SCRIPT VARIABLES
#
#
###############################################################
###############################################################
#
#
# Commented out
#
#
###############################################################
###############################################################

for value in filedata
	
	if (value[1]=="N_tau")	
	global N_tau = eval(value[2])
	end


	if(value[1]=="L") 
	global	L = eval(value[2])
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
	println("hello world")
	global file_name = String(eval(value[2]))
	end

end
println(file_name)
##############################################################################
##############################################################################






Na = L




####################################################
####################################################





# system size
N = L

# cutoff
#cutoff = 10^(-14)

# bond dimension
#maxdim=200

# type of spin used by mps and mpo
s = siteinds("S=1/2", N; conserve_qns= false)
#dt = 0.001*2*π
tau = tau*2*π

dt = tau/N_tau

#tau = tau*2*π
periodic = false
if N < 9
periodic = true
end
#total time parameters
ttotal = numkicks*tau
#hz = 0.1
    
gates = ITensor[]


# parameters kick strength
#phi_0 = 2*π*10^(-1)
#On site kick
kick= []

for site_and_phi in sites_and_phis 
	
	i = Int64(site_and_phi[1])
	phi = 2*π*site_and_phi[2]
	
		if phi > 10^(-14)
			OPKICK = cos(phi/2)*op("Id", s[i])+2*im*sin(phi/2)*op("Sx", s[i])
			push!(kick, OPKICK)
		end

end

###############################################################################
#
# Set up the hamiltonian and time evolution operator
#
################################################################################ 
	for j in 1:(N - 1)

		local i = Int64(bonds_and_Js[j][1])
		local ip = Int64(bonds_and_Js[j][2])

		J = bonds_and_Js[j][3]
		Jz = bonds_and_Js[j][4]
		
		local hz = sites_and_hzs[j][2]


		s1 = s[i]
		s2 = s[ip]
            
		hj = Jz* op("Sz", s1) * op("Sz", s2) + J* 1/2 * op("S+", s1) * op("S-", s2) +J * 1/2 * op("S-", s1) * op("S+", s2)- hz*op("Sz",s1)*op("Id",s2)

		Gj = exp(-im * dt/2 * hj)
		push!(gates, Gj)
       
	end


		J=bonds_and_Js[N][3]
		Jz = bonds_and_Js[N][4]

		hz = sites_and_hzs[N][2]
		s1 = s[N]
        	s2 = s[1]


	if periodic == true
		        
        	hj =Jz * op("Sz", s1) * op("Sz", s2) +J* 1/2 * op("S+", s1) * op("S-", s2) + J* 1/2 * op("S-", s1) * op("S+", s2) - hz*op("Sz",s1)*op("Id",s2)
    
       		Gj = exp(-im * dt/2 * hj)
        	push!(gates, Gj)
	

	else 
		hj =  - hz*op("Sz",s1)*op("Id",s2)
	 
       		Gj = exp(-im * dt/2 * hj)
        	push!(gates, Gj)
	
	end


    
	append!(gates,reverse(gates))


############################################################


############################################################


# Initialize fully polarized z state 
# first must make an array of string 
# up states

#############################################################

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

	global psi0 = apply(OP1,psi0;cutoff=cutoff1,maxdim =maxdim1)
	normalize(psi0)
	global psi0 = apply(OP2,psi0;cutoff=cutoff1,maxdim=maxdim1)
	normalize(psi0)
end


# make two MPSs one which is all X polarized state
# another which will be time evolved and one which will
# remain |ψ(0)> 
# here |ψ₁> ≡ exp(iσₓ ϕ/2) exp(-i H₀ τ) |ψ(0)> = |ψ(τ)>
# and  |ψ₂> = |ψ(0)>
 
psi1 = deepcopy(psi0)
psi2 = deepcopy(psi1)

println(" |⟨ψ|ψ⟩|² = ",abs2(inner(psi1,psi2)))

#########################################################
# initialize empty arrays for observables and correlators
#########################################################

time = []
Prob = []
Bonddim =[]
XXmat = []
Xmat = []
entropies = []


#########################################################
# Calculation of initial conditions for
# Echo time
#########################################################


# t₀ = 0
global t =0.0 



	# calculation of Echo at t₀ = 0
	Echo = abs2(inner(psi1,psi2))

	# time t₀ = 0
	push!(time,t)
	
	# Echo t₀ = 0
	push!(Prob,Echo)
	
	# max bond dimension at t₀ = 0
	println("t = ","maxbonddim = ", maxlinkdim(psi1))

	# SᵢˣSⱼˣ correlation matrix at t₀ = 0
	xxcorr= correlation_matrix(psi1,"Sx","Sx")
	push!(XXmat,xxcorr)

    # Sᵢˣ local observable at t₀= 0
    magz = expect(psi1,"Sx")
    push!(Xmat,magz)
	



############################################################
# Entropy Calculation at t = 0
# 
############################################################

	# initialize empty array of 
	# entropies will contain 
	# [[S₁(t₁),S₂(t₁),...,Sₙ(t₁)],[S₁(t₂),S₂(t₂),...,Sₙ(t₂)],
	#  [...],...,[S₁(tₙ),S₂(tₘ),...,Sₙ(tₙ)]],
	entropies = []
	
	
	entropy = []
	for σ in 1:N
		push!(entropy , entropy_von_neumann(psi1,σ))
	end 	
	
	push!(entropies, entropy)

############################################################
# End of Entropy 
# Calculation
############################################################
	








#########################################################
#
#
#########################################################

# Calculation for total time
# dt = τ / N_tau
# ttotal = numkicks τ
# ttotal = numkicks N_tau dt

for n in 0:N_tau*numkicks
	
	println(to)
	
	t≈ttotal && break

	# time evolution with H₀ without kick	
	global psi1 = apply(gates, psi1; cutoff = cutoff1,maxdim=maxdim1)
	normalize(psi1)
		
		# at t = n τ we kick the system. If statement to ensure
		# we are at an integral multiple of τ
		if n%N_tau == 0.0 && n > 0
			# kick
			for foot in kick
			psi1= apply(foot,psi1;cutoff=cutoff1,maxdim=maxdim1)
			normalize(psi1)
			end  
					# Echo at t = n τ
					local Echo = abs2(inner(psi1,psi2))
					
					# time t = nτ 
					push!(time,t)
					
					# Append array with Echo
					push!(Prob,Echo)
					
					# ξ is the max bond dim at t = n τ
					ξ =  maxlinkdim(psi1)
					push!(Bonddim, ξ)
					println("t = ","maxbonddim = ", ξ)

					# Sᵢˣ Sⱼˣ
					local xxcorr= correlation_matrix(psi1,"Sx","Sx")
					push!(XXmat,xxcorr)
         			       
					# Sᵢˣ local observable at t = n τ
					global magz = expect(psi1,"Sx")
        	        		push!(Xmat,magz)
	
				
				###############################################################	
				#  Beginning of Entropy 
				#  Calculation
				#  
				# 
				###############################################################

					local entropy = []
					for σ in 1:N
						push!(entropy , entropy_von_neumann(psi1,σ))
					end 	
	
					push!(entropies, entropy)

				###############################################################
				# End of Entropy
				# Calculation
				#
				#
				###############################################################


            	end
            
    	global t =dt+ t

end

time = time/(2*pi)



@save file_name time Prob Bonddim Xmat XXmat entropies









