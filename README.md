
<font face = "Times New Roman"> 
  
# Square Pulse Heisenberg TEBD
  
We implement time evolution for the staggered heisenberg model in a static Zeeman field with pulse transverse transverse field using matrix product states and Time Evolution Block Decimation.
 

$$H = H_0 + H_D$$

where
$$H_0 = \sum_{\langle i,j \rangle} (-1)^i \vec{S}_i \cdot \vec{S}_j +\sum_i h_z S_i^z \qquad H_D = {\rm sgn }\Big[{\rm sin \, \omega t}\Big]$$

To do this we have two hamiltonians

$$ H =  \begin{array}{ll}
      H_1 = H_0 +\sum_i \gamma_i S_i^x & {\rm for}\quad  0 \leq t\leq \tau/2 \\
       H_2 = H_0 -\sum_i \gamma_i S_i^x & {\rm for}\quad \tau/2 \leq t \leq \tau \\
   \end{array}
 $$

$\gamma_i$ is the strength of the transverse field on site $i$
  
We do time evolution with $H_1$ followed by time evolution with $H_2$.
  
$$F(\tau)  = e^{-i H_2 {\hspace{0.5cm}} \tau/2}  e^{-i H_1 {\hspace{0.5cm}} \tau/2}$$ 






