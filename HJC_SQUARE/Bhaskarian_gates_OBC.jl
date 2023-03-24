## Ronald Melendrez
# 1/23/2023
# SQUARE PULSE TEBD
using ITensors
using JLD2
using TimerOutputs

################################################
#
# GATES FUNCTION
#
################################################
include("swapFunction.jl")


function square_pulse_gate_generator(s,N,bonds_and_Js,sites_and_phis,sites_and_hzs,dt)


    ##########################################################
    # Initializing Gates
    ##########################################################


    # gates 1 for magnetic field along -X direction
    gates1 = ITensor[]

    # gates 2 for magnetic field along +X direction
    gates2 = ITensor[]


    ##########################################################
    # Gates 1 Operator
    ##########################################################
    for j in 1:N

            ############################
            # Indices and Couplings
            ############################

            # Nearest Neighbor Indices
            local i = Int64(bonds_and_Js[j][1])
            local ip = Int64(bonds_and_Js[j][2])

            # Heisenberg Couplings
            local J = Float64(bonds_and_Js[j][3])
            local Jz = Float64(bonds_and_Js[j][4])

            # Zeeman field couplings
            local hz = sites_and_hzs[j][2]

            # Pulse field coupling
            local hx = sites_and_phis[j][2]

            # site indices SᵢSⱼ
            s1 = s[i]
            s2 = s[ip]

            #################################################
            #      Hamiltonian Construction
            #################################################
        # OPEN BOUNDARY CONITIONS
        # Bulk Nearest Neighbor Condition
        if i+1 == ip || i-1 == ip

                # Ising couplings Jz SᵢᶻSⱼᶻ
                hj = Jz*op("Sz",s1)*op("Sz",s2)

                # XY couplings (J Sᵢ⁺Sⱼ⁻ + J Sᵢ⁻Sⱼ⁺)/2
                hj += J * (1/2) * op("S+",s1) * op("S-",s2) + J * (1/2) * op("S-",s1) * op("S+",s2)


                hj += -hz * (0.5 * op("Sz",s1) * op("Id",s2)+ 0.5 * op("Sz",s2) * op("Id",s1))



                # X pulse (same subtlety as the static Z field)
               hj += hx * op("Sx", s1) * op("Id", s2) 



        # OPEN BOUNDARY CONDITION
        else
                # bring 1 next to L-1 position so that it is next to L
                local U = swapSequence(s,ip,i)
                append!(gates1,U)
                s1 = s[i-1]
                s2 = s[i]

                # Ising couplings Jz SᵢᶻSⱼᶻ
                hj = Jz*op("Sz",s1)*op("Sz",s2)

                # XY couplings (J Sᵢ⁺Sⱼ⁻ + J Sᵢ⁻Sⱼ⁺)/2
                hj += J * (1/2) * op("S+",s1) * op("S-",s2) + J * (1/2) * op("S-",s1) * op("S+",s2)


                hj += -hz * (0.5 * op("Sz",s1) * op("Id",s2)+ 0.5 * op("Sz",s2) * op("Id",s1))



                # X pulse (same subtlety as the static Z field)
                hj += hx * op("Sx", s1) * op("Id", s2)


        end

        Gj = exp(-im * (dt/2) * hj)

        push!(gates1,Gj)
    end
        append!(gates1, reverse(gates1))



    ##########################################################
    # Gates 2 Operator
    ##########################################################

    
    for j in 1:N

            ############################
            # Indices and Couplings
            ############################

            # Nearest Neighbor Indices
            local i = Int64(bonds_and_Js[j][1])
            local ip = Int64(bonds_and_Js[j][2])

            # Heisenberg Couplings
            local J = Float64(bonds_and_Js[j][3])
            local Jz = Float64(bonds_and_Js[j][4])

            # Zeeman field couplings
            local hz = sites_and_hzs[j][2]

            # Pulse field coupling
            local hx = sites_and_phis[j][2]

            # site indices SᵢSⱼ
            s1 = s[i]
            s2 = s[ip]

            #################################################
            #      Hamiltonian Construction
            #################################################
        # OPEN BOUNDARY CONITIONS
        # Bulk Nearest Neighbor Condition if i and ip are on the left or on right of each other!
        if i+1 == ip || i-1 == ip

                # Ising couplings Jz SᵢᶻSⱼᶻ
                hj = Jz*op("Sz",s1)*op("Sz",s2)

                # XY couplings (J Sᵢ⁺Sⱼ⁻ + J Sᵢ⁻Sⱼ⁺)/2
                hj += J * (1/2) * op("S+",s1) * op("S-",s2) + J * (1/2) * op("S-",s1) * op("S+",s2)


                hj += -hz * (0.5 * op("Sz",s1) * op("Id",s2)+ 0.5 * op("Sz",s2) * op("Id",s1))



                # X pulse (same subtlety as the static Z field)
                hj += -hx  * op("Sx", s1) * op("Id", s2) 
            
else
                # bring 1 next to L-1 position so that it is next to L
                local U = swapSequence(s,ip,i)
                append!(gates2,U)
                s1 = s[i-1]
                s2 = s[i]

                # Ising couplings Jz SᵢᶻSⱼᶻ
                hj = Jz*op("Sz",s1)*op("Sz",s2)

                # XY couplings (J Sᵢ⁺Sⱼ⁻ + J Sᵢ⁻Sⱼ⁺)/2
                hj += J * (1/2) * op("S+",s1) * op("S-",s2) + J * (1/2) * op("S-",s1) * op("S+",s2)


                hj += -hz * (0.5 * op("Sz",s1) * op("Id",s2)+ 0.5 * op("Sz",s2) * op("Id",s1))



                # X pulse (same subtlety as the static Z field)
                hj += -hx * op("Sx", s1) * op("Id", s2) 


        end

        Gj = exp(-im * (dt/2) * hj)

        push!(gates2,Gj)
    end
        append!(gates2, reverse(gates2))
        return gates1, gates2

end
##########################################################################
#
# END OF GATES FUNCTION
#
#########################################################################


