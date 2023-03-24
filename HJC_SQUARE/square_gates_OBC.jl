## Ronald Melendrez
# 1/23/2023
# SQUARE PULSE TEBD
using ITensors
using JLD2
using TimerOutputs

################################################
# GATES FUNCTION
################################################
function square_pulse_gate_generator(s,N,bonds_and_Js,sites_and_phis,sites_and_hzs,dt)
    
    ##########################################################
    # Initializing Gates
    ##########################################################
    # gates 1 for magnetic field along -X direction
    gates1 = ITensor[]
    # gates 2 for magnetic field along +X direction
    gates2 = ITensor[]
    
    for case in 1:2
	    gammasign=1.0
	    if (case==2)
		gammasign=-1.0
	    end
	    for j in 1:(N-1)
	      ############################
	      # Indices and Couplings
	      ############################
	      # Nearest Neighbor Indices
	      local i = min(Int64(bonds_and_Js[j][1]),Int64(bonds_and_Js[j][2]))
	      local ip = max(Int64(bonds_and_Js[j][1]),Int64(bonds_and_Js[j][2]))

	      # Heisenberg Couplings
	      local J = Float64(bonds_and_Js[j][3])
	      local Jz = Float64(bonds_and_Js[j][4])

	      # Zeeman field couplings
	      local hzi = sites_and_hzs[i][2]
	      local hzip = sites_and_hzs[ip][2]

	      # Pulse field coupling
	      local gammaxi  = sites_and_phis[i][2]*gammasign
	      local gammaxip = sites_and_phis[ip][2]*gammasign

	      # site indices SᵢSⱼ
	      s1 = s[i]
	      s2 = s[ip]

	      #################################################
	      #      Hamiltonian Construction
	      #################################################
	      # OPEN BOUNDARY CONDITIONS
	      # Bulk Nearest Neighbor Condition
	      if i+1 == ip
			# Ising couplings Jz SᵢᶻSⱼᶻ
			hj = Jz*op("Sz",s1)*op("Sz",s2)

			# XY couplings (J Sᵢ⁺Sⱼ⁻ + J Sᵢ⁻Sⱼ⁺)/2
			hj += J * (1/2) * op("S+",s1) * op("S-",s2) + J * (1/2) * op("S-",s1) * op("S+",s2)

			if (i>1 && ip<N)
				hj += -hzi * 0.5 * op("Sz",s1) * op("Id",s2) - hzip* 0.5 * op("Id",s1) * op("Sz",s2)
				hj += gammaxi * 0.5 * op("Sx",s1) * op("Id",s2) + gammaxip* 0.5 * op("Id",s1) * op("Sx",s2)
			end

			if (i==1)
				hj += -hzi * 1.0 * op("Sz",s1) * op("Id",s2) - hzip* 0.5 * op("Id",s1) * op("Sz",s2)
				hj += gammaxi * 1.0 * op("Sx",s1) * op("Id",s2) + gammaxip* 0.5 * op("Id",s1) * op("Sx",s2)
			end

			if (ip==N)
				hj += -hzi * 0.5 * op("Sz",s1) * op("Id",s2) - hzip* 1.0 * op("Id",s1) * op("Sz",s2)
				hj += gammaxi * 0.5 * op("Sx",s1) * op("Id",s2) + gammaxip* 1.0 * op("Id",s1) * op("Sx",s2)
			end

	      end
	      Gj = exp(-im * (dt/2) * hj)
	      if (case==1)
			push!(gates1,Gj)
	      end
	      if (case==2)
			push!(gates2,Gj)
	      end
	    end
	    
	    if (case==1)    
		append!(gates1, reverse(gates1))
	    end 
	 
	    if (case==2)
		append!(gates2, reverse(gates2))
	    end
    end
    return gates1, gates2
end
##########################################################################
# END OF GATES FUNCTION
#########################################################################


