<<<<<<< HEAD
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

N = 10
s = siteinds("S=1/2", N)
psi = randomMPS(s, 4)
SvN = entropy_von_neumann(psi, b)
=======
"""I don't have thos code"""

"Prakash I believe in you. Lets make magnetism great again -Ron"
>>>>>>> a431a4a82f6c1774634061b9240c8afe28ce7e47
