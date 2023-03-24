# to run this script
# julia PLOTscript   filename                                                                                                                                                               System Size
# julia PLOT_TEBD.jl OPEN_Bhaskarian_pulse_pulse_40_site_1_site_8_output_with_effective_ham_double_S_experiment_1_site_1_kicked_L_8__hz_1_phi_15_tau_0.83442_maxdim_40_numkicks_10_Ntau_100 8

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
#using Plots
#using GR
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





# use data collected from file name to generate plot titles

# prints out the information letting the user
# know that such information has been extracted


# loads the information from the JDL2 file and stores into the following variables
@load file_name time Prob Bonddim Xmat XXmat entropies

#time = time*parse(Float64,tau)
# Some eigenstate overlaps are numerically zero.
# Because of this it is of value to record
# the number of non zero overlaps. This will
# provide a sanity check. HighE will be a list of
# non-zero H0 eigenstate probabilities
#time = time/tauvariable


println("time[1] = ", time[1],"Echo[1]=", Prob[1])
#pop!(Prob)
println("dim time = ", length(time),"dim Echo=", length(Prob))
#STOP
###########################################
# Echo plot  |<ψ(0)|ψ(t)>|^2
###########################################


#prakash_time = []
#prakash_Cos =[]
#x_prakash = 2\pi/100

#for i in 0:800
#	push!(prakash_time,x_prakash*i)
#	push!(prakash_Cos,(1+cos(x_prakash*i))/2)
#end




newtime = time
plot(2*pi* time/(2*pi/7.53),Prob)
#plot(prakash_time,prakash_Cos)
title(" Echo 1 D Staggered MG Spin-1/2\n")
ylabel(L"$|{\langle \psi(0)|\psi(t) \rangle}|^2$")
ylim(0,1.1)
savefig("Echo_PLOT"*fname*".png")
show()
close()

#timeprakash = []
#Echoprakash = []
#for (k,deltat) in enumerate(time)
#	push!(timeprakash,deltat)
#	push!(Echoprakash,Prob[k])
#	if k == 400
#		break
#	end
#
#nd

#plot(timeprakash,Echoprakash)
#how()






L =ARGS[2]
listoflegend = []
for ψ in 1:parse(Int64,L)
	push!(listoflegend,"site = "*string(ψ))
end
println(listoflegend)

#for t in time 
#	println("tn = ",t)
#end

println("S_10= ",entropies[10])
L_a=length(entropies[10])
chain = []
for i in 1:L_a
push!(chain,i)
end
x_L = length(chain)
y_L = length(time)
data = zeros(Float64,x_L,y_L)
xs = []
ys = []
zs = []
for x_i in 1:x_L
 for y_i in 1:y_L
	push!(xs,x_i)
	push!(ys,y_i)
	push!(zs,entropies[y_i][x_i])
end
end 



pcolormesh(entropies)
title(" Single Site Entropy Space Time Plot\n")
#ax[:locator_params](axis ="x", nbins=12)
xlabel("sites")
ylabel("time")
xlim(0,x_L+1)
clim(0,log(2))
colorbar()
legend()
#grid("on")
show()



plot(time/(2*pi/7.53),entropies)
title(" Single Site Entropy Spin-1/2\n")
axhline(log(2),color=:red)
ylim(0,1)
ylabel(L"$S_{VN}$")
legend(listoflegend)
savefig("SVN_"*fname*".png")
show()
close()

println("Svn site 1 1000 step = ",entropies[10][1])
println("Sx site 1 1000 step = ", Xmat[10][1])
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

plot(newtime/(2*pi/7.53),S1x)
legend(listoflegend)
end
title(" On-Site Spin Obs MG Spin-1/2\n")

ylabel(L"$\langle S_1^x \rangle$")
#ylabel(L"$|{\langle \psi(0)|\psi(t) \rangle}|^2$")
ylim(-0.5,0.5)
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
show()
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


#sara = Lby2-3
#for j in 1:parse(Int64,L)	
#SixSjx = []
#newtime = []
#for (i,element) in enumerate(XXmat)
#push!(newtime,time[parse(Int64,N_tau)*i]) 
#push!(SixSjx,real(element[sara,j]))
#end

###########################################
# Echo plot  |<ψ(0)|ψ(t)>|^2
###########################################


#subplot(311)
#axhline(1.0,color=:red)

#plot(newtime,SixSjx)
#legend(listoflegend)
#end
#title(" Corr MG Spin-1/2\n")

#ylabel(L"$\langle S_i^x S_j^x \rangle$"*" i = "*string(sara))
#savefig("S_"*string(sara)*"x_S_jx"*fname*".png")
#show()
#close()




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







#sara=Lby2-2
#for j in 1:parse(Int64,L)	
#SixSjx = []
#newtime = []
#for (i,element) in enumerate(XXmat)
#push!(newtime,time[parse(Int64,N_tau)*i]) 
#push!(SixSjx,real(element[sara,j]))
#end

#plot(newtime,SixSjx)
#legend(listoflegend)
#end
#title(" Corr MG Spin-1/2\n")
#ylabel(L"$\langle S_i^x S_j^x \rangle$"*" i = "*string(sara))
#savefig("S_"*string(sara)*"x_S_jx"*fname*".png")
#show()
#close()

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



