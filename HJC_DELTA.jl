using ITensors
using JLD2
using TimerOutputs
#using PyPlot
#using LaTeXStrings
using PrettyTables
using Printf

const to = TimerOutput()

include("single_site_entropy.jl")

####################################################
# function that reads in input data and stores it
function read_input_file(fname)
	data = []
	f = open(fname,"r")
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
# EXTRACTION OF SCRIPT DATA AND ASSIGNED TO SCRIPT VARIABLES
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
	println("bonds_and_Js = ",bonds_and_Js) 
	end
	
	if(value[1]=="maxdim") 
	global	maxdim1 = eval(value[2])
	end
	
	if(value[1]=="cutoff") 
	global	cutoff1 = eval(value[2])
	end
	
	if(value[1]=="sites_and_hzs") 
	global	sites_and_hzs = eval(value[2])
	println("sites_and_hzs = ",sites_and_hzs) 
	end

	if(value[1]=="sites_and_phis") 
	global	sites_and_phis = eval(value[2])
	println("sites and phis/2pi = ",sites_and_phis) 
	end
	
	if(value[1]=="tau")
	global	tau = eval(value[2])
	println("tau/2pi = ",tau) 
	end

	
	if(value[1]=="numkicks")
	global numkicks = eval(value[2])
	println("numkicks = ",numkicks) 
	end
	
	if(value[1]=="file_name")
	global file_name = String(eval(value[2]))
	end

end
println(file_name)
##############################################################################
##############################################################################

# system size
N = L

# type of spin used by mps and mpo
s = siteinds("S=1/2", N; conserve_qns= false)
tau = tau*2*π
dt = tau/N_tau
ttotal = numkicks*tau
periodic = false
gates = ITensor[]

# parameters kick strength
#On site kick
kick= []

for site_and_phi in sites_and_phis 
	i = Int64(site_and_phi[1])
	phi = 2*π*site_and_phi[2]
		if abs(phi) > 10^(-14)
			OPKICK = cos(phi/2)*op("Id", s[i])+2*im*sin(phi/2)*op("Sx", s[i])
			push!(kick, OPKICK)
		end

end

###############################################################################
# Set up the hamiltonian and time evolution operator
################################################################################ 
	for j in 1:(N - 1)  # works only for OBC

		local i = Int64(bonds_and_Js[j][1])
		local ip = Int64(bonds_and_Js[j][2])

		J = bonds_and_Js[j][3]
		Jz = bonds_and_Js[j][4]
		
		local hz = sites_and_hzs[j][2]

		s1 = s[i]
		s2 = s[ip]
            
		if (i==1 && ip==2)
			hj = Jz* op("Sz", s1) * op("Sz", s2) + J* 1/2 * op("S+", s1) * op("S-", s2) +J * 1/2 * op("S-", s1) * op("S+", s2)- hz*op("Sz",s1)*op("Id",s2) - (0.5*hz)*op("Id",s1)*op("Sz",s2)
		elseif (i==N-1 && ip==N)
			hj = Jz* op("Sz", s1) * op("Sz", s2) + J* 1/2 * op("S+", s1) * op("S-", s2) +J * 1/2 * op("S-", s1) * op("S+", s2)- (0.5*hz)*op("Sz",s1)*op("Id",s2) - hz*op("Id",s1)*op("Sz",s2)
		else
			hj = Jz* op("Sz", s1) * op("Sz", s2) + J* 1/2 * op("S+", s1) * op("S-", s2) +J * 1/2 * op("S-", s1) * op("S+", s2)- (0.5*hz)*op("Sz",s1)*op("Id",s2) - (0.5*hz)*op("Id",s1)*op("Sz",s2)
		end
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
#  Polarize spins along initial directions
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
for n in 0:N_tau*numkicks
	global t=n*dt
	t≈ttotal && break
	##########################################################
	# Do measurement
	##########################################################
	# at t = integer τ we kick the system. If statement to ensure
	# we are at an integral multiple of τ
	if n%N_tau == 0  ## then do measurement 
		# Echo at t = integer τ
		local Echo = abs2(inner(psi1,psi2))
		# ξ is the max bond dim at t = integer τ
		ξ =  maxlinkdim(psi1)
		#push!(Bonddim, ξ)
		@printf("Kick = %4d t = %+5.15f Echo = %+5.15f  maxbonddim = %4d \n",n/N_tau,t,Echo,ξ)
		# Sᵢˣ Sⱼˣ
		#local xxcorr= correlation_matrix(psi1,"Sx","Sx")
		# Sᵢˣ local observable at t = n τ
		local magx = expect(psi1,"Sx")
		local magy = expect(psi1,"Sy")
		local magz = expect(psi1,"Sz")
		###############################################################	
		#  Beginning of Entropy Calculation
		###############################################################
		println("Site        <Sx>                   <Sy>                <Sz>                   Svn")
		for σ in 1:N
			@printf("%3d  %+5.15f  %+5.15f  %+5.15f  %+5.15f \n", σ,magx[σ],magy[σ],magz[σ],entropy_von_neumann(psi1,σ))
		end
		println(" ") 	
		###############################################################
		# End of Entropy Calculation
		###############################################################
	end # measurement

	##########################################################
	# Time evolution
	###########################################################
	# time evolution with H₀ without kick exp(-i Ho dt)	
	global psi1 = apply(gates, psi1; cutoff = cutoff1,maxdim=maxdim1)
	normalize(psi1)		
	# at t = integer τ we kick the system. If statement to ensure
	# we are at an integral multiple of τ
	if (n+1)%N_tau == 0
			# kick
			for foot in kick
				psi1= apply(foot,psi1;cutoff=cutoff1,maxdim=maxdim1)
				normalize(psi1)
			end  
	end
end
