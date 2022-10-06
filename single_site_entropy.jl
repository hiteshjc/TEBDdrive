using ITensors

function entropy_von_neumann(psi::MPS, b::Int)
  s = siteinds(psi)  
  orthogonalize!(psi, b)
  psidag = dag(psi)
  rho = psi[b]*prime(psidag[b],"Site")
  S,_ = diagonal(rho) """ replace with correct julia command """
  SvN = 0.0
  for n in 1:dim(S, 1)
    p = S[n,n]
    SvN -= p * log(p)
  end
  return SvN
end

N = 4
s = siteinds("S=1/2", N)
psi = randomMPS(s, 4)
b = 2
SvN = entropy_von_neumann(psi, b)