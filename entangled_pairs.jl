using ITensors
cutoff = 10^(-14)
N =3
entangled_pairs = [[1,2]]
s = siteinds("S=1/2", N; conserve_qns= false)
maxdim = 200

include("single_site_entropy.jl")
include("two_site_RDM_n_SvN.jl")



function singletgenerator(N,entangled_pairs)
states = []
	for n in 1:N
        	 push!(states, "Up")
	end


	global psi0 = MPS(s,states)
	
	for entangled_pair in entangled_pairs
        	 i = entangled_pair[1]
        	 j = entangled_pair[2]

         	OPsinglet =op("Id", s[i])*op("S-",s[j])-op("Id", s[j])*op("S-",s[i])

         	global psi0 = apply(OPsinglet,psi0;cutoff,maxdim)
         	psi0 = normalize(psi0)
	end


return psi0

end

#################################################################
# ONE SITE ENTANGLEMENT MEASUREMENTS
#
#################################################################


psi0 = singletgenerator(N,entangled_pairs)


println(" (↑₁↓₂ - ↓₁↑₂)/√2     S₁ = ", entropy_von_neumann(psi0, 1))
println(" (↑₁↓₂ - ↓₁↑₂)/√2     S₂ = ", entropy_von_neumann(psi0, 2))
states = []
	for n in 1:N
        	 push!(states, "Up")
	end


	global psi1= MPS(s,states)
	

println(" ↑₁↑₂     S₁ = ", entropy_von_neumann(psi1, 1))


########################################################
# TWO SITE ENTANGLEMENT MEASUREMENTS
#
########################################################


println("(↑₁↓₂↑₃-↓₁↑₂↑₃) / √2       S₁₃ = ",two_site_entropy_von_neumann(psi0,1,3)," should be S₁₃ = ln(1) =0 because ")
println("(↑₁↓₂↑₃-↓₁↑₂↑₃) / √2       S₁₂ = ",two_site_entropy_von_neumann(psi0,1,2)," should be =  ln 2² = 1.38629436112")

