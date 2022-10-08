using LinearAlgebra
using ITensors

#input mush have b < c
function twoSite_entropy_von_neumann(psi::MPS, b::Int, c::Int)
  if b > c 
	throw("b = $b must be less than c = $c")
  end
  
  s = siteinds(psi)  
  orthogonalize!(psi, b)
  psidag = dag(psi)
  prime(psidag,"Link")

  #index linking b to b-1
  lb = commonind(psi[b],psi[b-1])
  rho = prime(psi[b], lb)*prime(psidag[b],"Site")
  for n in b+1:c-1
	rho *= psi[n]
	rho *= psidag[n]
  end

  #index linking c to c+1
  lc = commonind(psi[c],psi[c+1])
  rho *= prime(psi[c],lc)
  rho *= prime(psidag[c],"Site")
  
  S,_ = eigen(rho)
  SvN = 0.0
  for n in 1:dim(S, 1)
     p = real(S[n,n])
	if p > 10^-10
        SvN -= p * log(p)
	end
  end
  return SvN
  end 
