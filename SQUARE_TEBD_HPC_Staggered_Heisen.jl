using ITensors
using JLD2
using TimerOutputs
#using PyPlot
#using LaTeXStrings
using PrettyTables

const to = TimerOutput()

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

##
#L = 2
#numkicks = 1000
# kick period τ
##
#tau = 2*π*0.01
for value in filedata
	
	if (value[1]=="N_tau")	
	global N_tau = eval(value[2])
	end


	if(value[1]=="L") 
	global	L = eval(value[2])
	end
##	
#hz= 0.1	
#sites_and_initial_product_states = [[1,[1/sqrt(2),0.0],[1/sqrt(2),pi]],[2,[1/sqrt(2),0.0],[1/sqrt(2),pi]]]

##
#hz = 0.1
#sites_and_hzs=[[1,hz],[2,hz]]
	if(value[1]=="sites_and_initial_product_states") 
	global	sites_and_initial_product_states = eval(value[2])
	end

##
#bonds_and_Js = [[1,2,1.0,1.0],[1,2,-1.0,-1.0]]

# parameters kick strength

##
#phi_0 = 2*π*10^(-1)
#sites_and_phis=[[1,phi_0],[2,0.0],[3,0.0],[4,0.0],[5,0.0],[6,0.0],[7,0.0],[8,0.0],[9,0.0],[10,0.0]]
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






#Na =2

Na = L




####################################################
####################################################





# system size
N = L

# cutoff
cutoff1 = 10^(-14)

# bond dimension
#maxdim=200

# type of spin used by mps and mpo
s = siteinds("S=1/2", N; conserve_qns= false)
#dt = 0.001*2*π
tau = tau*2*π

# tau = 2π *0.1
dt = tau/N_tau

#tau = tau*2*π
periodic = false
if N < 9
periodic = true
end
#total time parameters
ttotal = numkicks*tau
#hz = 0.1

#gates 1
gates1 = ITensor[]

#gates 2
gates2 = ITensor[]



# parameters kick strength
#phi_0 = 2*π*10^(-1)
#On site kick
kick= []



##################################################################################
# KICK OPERATOR CONSTRUCTION
#
##################################################################################

for site_and_phi in sites_and_phis 
	
	i = Int64(site_and_phi[1])
	phi = 2*π*site_and_phi[2]
	
		if phi > 10^(-14)
			OPKICK = cos(phi/2)*op("Id", s[i])+2*im*sin(phi/2)*op("Sx", s[i])
			push!(kick, OPKICK)
		end

end


###############################################################################
# TIME EVOLUTION OPERATOR CONSTRUCTION
#
###############################################################################

###############################################################################
# GATES1
# Set up the hamiltonian and time evolution operator
# Magnetic Field in the negative x direciton
################################################################################ 
	for j in 1:(N - 1)

		local i = Int64(bonds_and_Js[j][1])
		local ip = Int64(bonds_and_Js[j][2])
		
		local iphi = Int64(sites_and_phis[j][1])		
		local hx = sites_and_phis[j][2]

		J = bonds_and_Js[j][3]
		Jz = bonds_and_Js[j][4]
		
		local hz = sites_and_hzs[j][2]


		s1 = s[i]
		s2 = s[ip]
            
		hj = Jz* op("Sz", s1) * op("Sz", s2) + J* 1/2 * op("S+", s1) * op("S-", s2) +J * 1/2 * op("S-", s1) * op("S+", s2)- hz*op("Sz",s1)*op("Id",s2) + hx*op("Sx",s1)*op("Id",s2)

		Gj = exp(-im * (dt/2) * hj)
		push!(gates1, Gj)
       
	end


		J=bonds_and_Js[N][3]
		Jz = bonds_and_Js[N][4]

		hz = sites_and_hzs[N][2]
		s1 = s[N]
        	s2 = s[1]
		hx = sites_and_phis[N][2]

	if periodic == true
			        
        	hj =Jz * op("Sz", s1) * op("Sz", s2) +J* 1/2 * op("S+", s1) * op("S-", s2) + J* 1/2 * op("S-", s1) * op("S+", s2) - hz*op("Sz",s1)*op("Id",s2) + hx*op("Sx",s1)*op("Id",s2)
    
       		Gj = exp(-im * (dt/2) * hj)
        	push!(gates1, Gj)
	

	else 
		hj =  - hz*op("Sz",s1)*op("Id",s2) + hx*op("Sx",s1)*op("Id",s2)
	 
       		Gj = exp(-im * (dt/2) * hj)
        	push!(gates1, Gj)
	
	end


    
	append!(gates1,reverse(gates1))

###############################################################################
# GATES2
# Set up the hamiltonian and time evolution operator
# Magnetic Field in the positive X direciton
################################################################################ 
	for j in 1:(N - 1)

		local i = Int64(bonds_and_Js[j][1])
		local ip = Int64(bonds_and_Js[j][2])
		
		local iphi = Int64(sites_and_phis[j][1])		
		local hx = sites_and_phis[j][2]

		J = bonds_and_Js[j][3]
		Jz = bonds_and_Js[j][4]
		
		local hz = sites_and_hzs[j][2]


		s1 = s[i]
		s2 = s[ip]
            
		hj = Jz* op("Sz", s1) * op("Sz", s2) + J* 1/2 * op("S+", s1) * op("S-", s2) +J * 1/2 * op("S-", s1) * op("S+", s2)- hz*op("Sz",s1)*op("Id",s2) - hx*op("Sx",s1)*op("Id",s2)

		Gj = exp(-im * (dt/2) * hj)
		push!(gates2, Gj)
       
	end


		J=bonds_and_Js[N][3]
		Jz = bonds_and_Js[N][4]

		hz = sites_and_hzs[N][2]
		s1 = s[N]
        	s2 = s[1]
		hx = sites_and_phis[N][2]


	if periodic == true
		        
        	hj =Jz * op("Sz", s1) * op("Sz", s2) +J* 1/2 * op("S+", s1) * op("S-", s2) + J* 1/2 * op("S-", s1) * op("S+", s2) - hz*op("Sz",s1)*op("Id",s2) - hx*op("Sx",s1)*op("Id",s2)
    
       		Gj = exp(-im * (dt/2) * hj)
        	push!(gates2, Gj)
	

	else 
		hj =  - hz*op("Sz",s1)*op("Id",s2) - hx*op("Sx",s1)*op("Id",s2)
	 
       		Gj = exp(-im * (dt/2) * hj)
        	push!(gates2, Gj)
	
	end


    
	append!(gates2,reverse(gates2))


############################################################



##############################################################################
# PREPARING INITIAL STATE
#
##############################################################################

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

println("Hello world = ",abs2(inner(psi1,psi2)))

time = []
Prob = []
Bonddim =[]
XXmat = []
Xmat = []
#for t in 0.0:dt:ttotal
global t =0.0 

# dt = τ / N_tau
# ttotal = numkicks τ
# ttotal = numkicks N_tau dt
n_tau = 1
	Echo = abs2(inner(psi1,psi2))
	println("|<ϕ|U|ϕ>|^2= ",Echo)
	push!(time,t)
	push!(Prob,Echo)
        #ξ=linkdims(psi1)
 	#println("bonddim = ",ξ )    
	push!(Bonddim, maxlinkdim(psi1))
	println("t = ","maxbonddim = ", maxlinkdim(psi1))

	xxcorr= correlation_matrix(psi1,"Sx","Sx")
	push!(XXmat,xxcorr)

	magz = expect(psi1,"Sx")
	push!(Xmat,magz)

switch = true
for n in 1:N_tau*numkicks
	

	println(to)

    	global t =dt+ t
    	t≈ttotal && break
	
	
                if n_tau <= N_tau/2
			global psi1 = apply(gates1, psi1; cutoff1,maxdim1)
			psi1=normalize(psi1)
	               end


                if n_tau > N_tau/2
 			global psi1 = apply(gates2, psi1; cutoff1,maxdim1)
			psi1 = normalize(psi1)
                     
                        if n_tau == N_tau
                         	push!(time,t)
				            xxcorr= correlation_matrix(psi1,"Sx","Sx")
				            push!(XXmat,xxcorr)
				            #println("<SᵢˣSⱼˣ> tables")
         			       
				            #pretty_table(xxcorr)
                		        #println(" ")
                

				            magz = expect(psi1,"Sx")
                		    push!(Xmat,magz)
				
			                local Echo = abs2(inner(psi1,psi2))
			                push!(Prob,Echo)
        		            #ξ=linkdims(psi1)
			                push!(Bonddim, maxlinkdim(psi1))
			                println("t = ","maxbonddim = ", maxlinkdim(psi1))


				            n_tau = 1
                            global switch = false
                        end

                end



                if switch ==  true
                        global  n_tau = 1 + n_tau
                end

                switch = true 


    end

time = time/(2*pi)


#file_name = "data_1"

for ele in Bonddim

println("Bonddim = ", ele)

end



@save file_name time Prob Bonddim Xmat XXmat










#plot(time,Prob)
#title("TEBD |X> state time evolution Staggered Heis + Zeeman hz = 0.1 " )
#ylim(0.0,1.1)
#xlabel(L"$2\pi/J$")
#ylabel(L"$|\langle \psi(0) | \psi(t) \rangle|^2$")
#show()
#x = []
#for i in 1:length(Bonddim)
#push!(x, Bonddim[i][1])
#end

#plot(time,x)
#plot(time,y)
#plot(time,z)
#show()



