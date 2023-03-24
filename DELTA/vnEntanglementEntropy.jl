using ITensors

function entropy_von_neumann(psi::MPS, b::Int)
  s = siteinds(psi)  
  orthogonalize!(psi, b)
  _,S = svd(psi[b], (linkind(psi, b-1), s[b]))
  SvN = 0.0
  for n in 1:dim(S, 1)
    p = S[n,n]^2
    SvN -= p * log(p)
  end
  return SvN
end

N = 4
s = siteinds("S=1/2", N)
psi = randomMPS(s, 4)
b = 2
SvN = entropy_von_neumann(psi, b)
@show SvN
