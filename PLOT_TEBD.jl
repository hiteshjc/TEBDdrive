# Annotated 22 Jan 2022
# 
# H0:    Time indenpendent Heisenberg
# H1(t): Periodic Tranverse Field
#
# H(t) =  H0 + H1(t)
#
# Scar Defects Project Time Evolution Plotting Script
# 1. Its argument is a JDL2 file type
# 2. Plots the overlap of |ψ(0)> onto the eigenstates of H0 
# 3. Plots Loschmidt echo of whole hamiltonian H(t)
# 4. Plots the overlap of |ψ(0)> onto eigenstates of F
#
##########################################################################


ENV["gkswstype"] = "gksqt"

using LaTeXStrings
using JLD2
using ArgParse
using PyPlot
global legends = []


##############################################
# Code Takes in name of file to be analyzed
##############################################

# input from terminal
fname = ARGS[1]
file_name =fname
# prints fname 
println(fname)

# filename contains information about the run
# here we extract such information
namepieces = split(fname,"_")
for (i,namepiece) in enumerate(namepieces)
         if namepiece =="L"
         global  L = string(namepieces[i+1])
         end
         if namepiece =="phi"
         global  phi = string(namepieces[i+1])
         end
         if namepiece == "tau"
         global        tau = string(namepieces[i+1])
         end
         if namepiece == "experiment"
         global      experiment = string(namepieces[i+1])
         end

         if namepiece == "maxdim"
         global      maxdim = string(namepieces[i+1])
         end

        if namepiece == "numkicks"
         global      numkicks = string(namepieces[i+1])
         end

         if namepiece == "Ntau"
         global      N_tau = string(namepieces[i+1])
         end

	if namepiece =="TEBD"
        global protocol = string(namepieces[i+1])
        end




end

# use data collected from file name to generate plot titles
plottitle1 = "experiment "*string(experiment)*"  "*L"$L = $"*string(L)
plottitle2= L"$\phi =$"*phi*"  "*L"$\tau =$"*tau*"  "*L"$h_z = 0.1$"*"  "*L"$maxdim = $"*maxdim*"  "*L"$N_\tau = $"*N_tau

# prints out the information letting the user
# know that such information has been extracted
println(plottitle1, "\n ",plottitle2)


# loads the information from the JDL2 file and stores into the following variables
@load file_name time Prob Bonddim Xmat XXmat entropies

#time = time*parse(Float64,tau)
# Some eigenstate overlaps are numerically zero.
# Because of this it is of value to record
# the number of non zero overlaps. This will
# provide a sanity check. HighE will be a list of
# non-zero H0 eigenstate probabilities
highE = []
tauvariable=parse(Float64, tau) 
#time = time/tauvariable


println("time[1] = ", time[1],"Echo[1]=", Prob[1])
#pop!(Prob)
println("dim time = ", length(time),"dim Echo=", length(Prob))
#STOP
###########################################
# Echo plot  |<ψ(0)|ψ(t)>|^2
###########################################
newtime = time
plot(time,Prob)
title(protocol*" Echo 1 D Chain Spin-1/2\n"*string(plottitle1)*"\n"*string(plottitle2))
ylabel(L"$|{\langle \psi(0)|\psi(t) \rangle}|^2$")
#ylim(-4,1.1)
savefig("Echo_PLOT"*fname*".png")
#show()
close()

listoflegend = []
for ψ in 1:parse(Int64,L)
	push!(listoflegend,"site = "*string(ψ))
end
println(listoflegend)
plot(time,entropies)
title(protocol *" Single Site Entropy Spin-1/2\n"*string(plottitle1)*"\n"*string(plottitle2))
axhline(log(2),color=:red)
ylim(0,1)
ylabel(L"$S_{VN}$")
legend(listoflegend)
savefig("SVN_"*fname*".png")
show()
close()

listoflegend = []
for k in 1:parse(Int64,L)
	push!(listoflegend, "j="*string(k))
end

legend(listoflegend)
global Lby2 = Int64(parse(Int64,L)/2)
#println("this is tau ", tau, " ", typeof(tau))
println("Xmat[1]=",Xmat[1])



###########################################
# On-Site Observable  <ψ(t)|Sᵢˣ|ψ(t)>
###########################################


#
#

for j in 1:parse(Int64,L)	
S1x = []
#newtime = []
for (i,element) in enumerate(Xmat)
#push!(newtime,time[parse(Int64,N_tau)*i]) 
push!(S1x,element[j])
end

#subplot(311)
#axhline(1.0,color=:red)

plot(newtime,S1x)
legend(listoflegend)
end
title(protocol*" On-Site Spin Obs 1 D Chain Spin-1/2\n"*string(plottitle1)*"\n"*string(plottitle2))

ylabel(L"$\langle S_1^x \rangle$")
#ylabel(L"$|{\langle \psi(0)|\psi(t) \rangle}|^2$")
#ylim(-4,1.1)
#subplot(312)	
#axhline(log(2),color=:red)
#plot(time,real.(entropies_1))
#ylabel("Entropy for Site 1")
#ylim(-0.1,1.1)
#axhline(log(2),color=:red)
#subplot(313)	
#plot(time,real.(entropies_2))
#axhline(log(2),color=:red)
#ylabel("Entropy for Site 5")
#xlabel("Number of kicks "*L"$\times 10^{3}$")
#ylim(-0.1,1.1)
savefig("S_ix_PLOT"*fname*".png")
#show()
close()
a = XXmat[1]
println("XXmat[1,1]= ", a[1,1])

###########################################
# Two Body Correlator  <ψ(t)|SᵢˣSⱼˣ|ψ(t)>  
#     i = L/2 -3  
#     j = [1,2,3,4,5,...,L]
#     
#     Suppose L = 10  kicked site will be site 5
#     We want to measure how correlated the spin 
#     3 lattice constants aways is to the 
#     other spins in the system.
#
#     . . . . . . . . . . 
#       ↑     ↑           
#       i     K           
###########################################


sara = Lby2-3
for j in 1:parse(Int64,L)	
SixSjx = []
#newtime = []
for (i,element) in enumerate(XXmat)
#push!(newtime,time[parse(Int64,N_tau)*i]) 
push!(SixSjx,real(element[sara,j]))
end

###########################################
# Echo plot  |<ψ(0)|ψ(t)>|^2
###########################################


#subplot(311)
#axhline(1.0,color=:red)

plot(newtime,SixSjx)
legend(listoflegend)
end
title(protocol*" Corr 1 D Chain Spin-1/2\n"*string(plottitle1)*"\n"*string(plottitle2))

ylabel(L"$\langle S_i^x S_j^x \rangle$"*" i = "*string(sara))
savefig("S_"*string(sara)*"x_S_jx"*fname*".png")
#show()
close()




###########################################
# Two Body Correlator  <ψ(t)|SᵢˣSⱼˣ|ψ(t)>  
#     i = L/2 -2  
#     j = [1,2,3,4,5,...,L]
#     
#     Suppose L = 10  kicked site will be site 5
#     We want to measure how correlated the spin 
#     2 lattice constants aways is to the 
#     other spins in the system.
#
#     . . . . . . . . . . 
#         ↑   ↑           
#         i   K           
###########################################







sara=Lby2-2
for j in 1:parse(Int64,L)	
SixSjx = []
#newtime = []
for (i,element) in enumerate(XXmat)
#push!(newtime,time[parse(Int64,N_tau)*i]) 
push!(SixSjx,real(element[sara,j]))
end

plot(newtime,SixSjx)
legend(listoflegend)
end
title(protocol*" Corr 1 D Chain Spin-1/2\n"*string(plottitle1)*"\n"*string(plottitle2))
ylabel(L"$\langle S_i^x S_j^x \rangle$"*" i = "*string(sara))
savefig("S_"*string(sara)*"x_S_jx"*fname*".png")
#/show()
close()

###########################################
# Two Body Correlator  <ψ(t)|SᵢˣSⱼˣ|ψ(t)>  
#     i = L/2 - 1  
#     j = [1,2,3,4,5,...,L]
#     
#     Suppose L = 10  kicked site will be site 5
#     We want to measure how correlated the spin 
#     1 lattice constant aways is to the 
#     other spins in the system.
#
#     . . . . . . . . . . 
#           ↑ ↑           
#           i K           
###########################################



sara=Lby2-1
for j in 1:parse(Int64,L)	
SixSjx = []
#newtime = []
for (i,element) in enumerate(XXmat)
#push!(newtime,time[parse(Int64,N_tau)*i]) 
push!(SixSjx,real(element[sara,j]))
end


plot(newtime,SixSjx)
legend(listoflegend)
end
title(protocol*" Corr 1 D Chain Spin-1/2\n"*string(plottitle1)*"\n"*string(plottitle2))
ylabel(L"$\langle S_i^x S_j^x \rangle$"*" i = "*string(sara))
savefig("S_"*string(sara)*"x_S_jx"*fname*".png")
#show()
close()

###########################################
# Two Body Correlator  <ψ(t)|SᵢˣSⱼˣ|ψ(t)>  
#     i = L/2 - 1  
#     j = [1,2,3,4,5,...,L]
#     
#     Suppose L = 10  kicked site will be site 5
#     We want to measure how correlated the spin 
#     1 lattice constant aways is to the 
#     other spins in the system.
#
#     . . . . . . . . . . 
#             ↑ ↑          
#             K i          
###########################################


sara=Lby2+1
for j in 1:parse(Int64,L)	
SixSjx = []
#newtime = []
for (i,element) in enumerate(XXmat)
#push!(newtime,time[parse(Int64,N_tau)*i]) 
push!(SixSjx,real(element[sara,j]))
end

plot(newtime,SixSjx)
legend(listoflegend)
end
title(protocol*" Corr 1 D Chain Spin-1/2\n"*string(plottitle1)*"\n"*string(plottitle2))
ylabel(L"$\langle S_i^x S_j^x \rangle$"*" i = "*string(sara))
savefig("S_"*string(sara)*"x_S_jx"*fname*".png")
#show()
close()
###########################################
# Two Body Correlator  <ψ(t)|SᵢˣSⱼˣ|ψ(t)>  
#     i = L/2 + 2  
#     j = [1,2,3,4,5,...,L]
#     
#     Suppose L = 10  kicked site will be site 5
#     We want to measure how correlated the spin 
#     2 lattice constants aways is to the 
#     other spins in the system.
#
#     . . . . . . . . . . 
#             ↑   ↑        
#             K   i         
###########################################




sara=Lby2+2
for j in 1:parse(Int64,L)	
SixSjx = []
#newtime = []
for (i,element) in enumerate(XXmat)
#push!(newtime,time[parse(Int64,N_tau)*i]) 
push!(SixSjx,real(element[sara,j]))
end

plot(newtime,SixSjx)
legend(listoflegend)
end
title(protocol*" Corr 1 D Chain Spin-1/2\n"*string(plottitle1)*"\n"*string(plottitle2))
ylabel(L"$\langle S_i^x S_j^x \rangle$"*" i = "*string(sara))
savefig("S_"*string(sara)*"x_S_jx"*fname*".png")
#show()
close()


sara=Lby2+3
for j in 1:parse(Int64,L)	
SixSjx = []
#newtime = []
for (i,element) in enumerate(XXmat)
#push!(newtime,time[parse(Int64,N_tau)*i]) 
push!(SixSjx,real(element[sara,j]))
end

plot(newtime,SixSjx)
legend(listoflegend)
end
title(protocol*" Corr 1 D Chain Spin-1/2\n"*string(plottitle1)*"\n"*string(plottitle2))
ylabel(L"$\langle S_i^x S_j^x \rangle$"*" i = "*string(sara))
savefig("S_"*string(sara)*"x_S_jx"*fname*".png")
#show()
close()












