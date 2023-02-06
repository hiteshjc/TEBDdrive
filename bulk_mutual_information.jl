# Feb 6 2023
# two site reduced density matrix 
# and mutual information calculation
#
using ITensors
function bult_mutual_info(pairs,psi,N)
    s = siteinds(psi)
    i = pairs[1]
    j = pairs[2]
    psi=orthogonalize(psi, i)
    psidag = dag(psi)
    prime!(psidag,"Link")
    


    lb = commonind(psi[i],psi[i-1])
    rho = prime(psi[i],lb)*prime(psidag[i],"Site")
    for k in i+1:j-1
        rho *= psi[k]
        rho *= psidag[k]
    end

    lc = commonind(psi[j],psi[j+1])

    rho *= prime(psi[j],lc)
    rho*= prime(psidag[j],"Site")

    
    rho *= psi[j+1]
    rho *= psidag[j+1]

    S,_ = eigen(rho)
    SvN = 0.0
  
     for n in 1:dim(S, 1)
        p = real(S[n,n])
            if abs(p) > 10^(-14)
                SvN -= round(p * log(p),digits=10)
            end
     end
    
      return SvN
end       
