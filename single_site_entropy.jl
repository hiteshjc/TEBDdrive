using LinearAlgebra
using ITensors

function entropy_von_neumann(psi::MPS, b::Int)
  s = siteinds(psi)  
  orthogonalize!(psi, b)
  psidag = dag(psi)
  rho = psi[b]*prime(psidag[b],"Site")
  S,_ = eigen(rho)
  SvN = 0.0
  for n in 1:dim(S, 1)
    p = real(S[n,n])
    SvN -= p * log(p)
  end
  return SvN
end

