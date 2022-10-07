using ITensors
cutoff = 10^(-14)
N =2
entangled_pairs = [[1,2]]
s = siteinds("S=1/2", N; conserve_qns= false)
maxdim = 200

include("single_site_entropy.jl")




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


psi0 = singletgenerator(N,entangled_pairs)
site = 1 
println(" (↑↓ - ↓↑)/√2     S = ", entropy_von_neumann(psi0, site))

states = []
	for n in 1:N
        	 push!(states, "Up")
	end


	global psi1= MPS(s,states)
	

println(" ↑↑     S = ", entropy_von_neumann(psi1, site))
