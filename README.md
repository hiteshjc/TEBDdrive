
<font face = "Times New Roman"> 

# Delta Kick and Pulse Field Heisenberg TEBD
  
We implement time evolution for the staggered heisenberg model in a static Zeeman field with both a delta kick and pulse transverse field protocols using matrix product states and Time Evolution Block Decimation.
 

$$H = H_0 + H_D$$

## Delta kick protocol
$$H_0 = \sum_{\langle i,j \rangle} (-1)^i \vec{S}_i \cdot \vec{S}_j +\sum_i h_z S_i^z \qquad H_D = \sum_n \sum_i \delta(t - n \tau) \phi^x_i S_i^x$$

$$F(\tau)  = \Big[\overset{N}{\underset{j=1}{\Pi}} e^{-i \frac{\phi_i}{2}  \sigma_j^x}\Big]  e^{-i H_0 {\hspace{0.5cm}} \tau/2}$$ 
## Pulse protocol
where
$$H_0 = \sum_{\langle i,j \rangle} (-1)^i \vec{S}_i \cdot \vec{S}_j +\sum_i h_z S_i^z \qquad H_D = {\rm sgn }\Big[{\rm sin \, \omega t}\Big]\sum_i \gamma_i S_i^x$$

To do this we have two hamiltonians

$$ H = \Bigg\lbrace \begin{array}{ll}
      H_1 = H_0 +\sum_i \gamma_i S_i^x & {\rm for}\quad  0 \leq t\leq \tau/2 \\
       H_2 = H_0 -\sum_i \gamma_i S_i^x & {\rm for}\quad \tau/2 \leq t \leq \tau \\
   \end{array}
 $$

$\gamma_i$ is the strength of the transverse field on site $i$
  
We do time evolution with $H_1$ followed by time evolution with $H_2$.
  
$$F(\tau)  = e^{-i H_2 {\hspace{0.5cm}} \tau/2}  e^{-i H_1 {\hspace{0.5cm}} \tau/2}$$ 

### Entanglement entropy and Mutual Information

$$|{\rm singlet} \rangle =\frac{1}{\sqrt{2}}(|\uparrow_1 \downarrow_2\rangle - |\downarrow_1 \uparrow_2 \rangle)$$

$$\rho = | \psi \rangle \langle \psi | = \frac{1}{\sqrt{2}}(|\uparrow_1 \downarrow_2\rangle - |\downarrow_1 \uparrow_2 \rangle) \cdot \frac{1}{\sqrt{2}}(\langle \uparrow_1 \downarrow_2| - \langle \downarrow_1 \uparrow_2 |)$$

$$\rho = | \psi \rangle \langle \psi | = \frac{1}{2}( |\uparrow_1 \downarrow_2 \rangle \langle \uparrow_1 \downarrow_2| - |\downarrow_1 \uparrow_2 \rangle \langle \uparrow_1 \downarrow_2| -|\uparrow_1 \downarrow_2 \rangle \langle \downarrow_1 \uparrow_2 | + |\downarrow_1 \uparrow_2 \rangle \langle \downarrow_1 \uparrow_2 |) $$

$$ \rho_1 = {\rm Tr_{2}}(\rho)  = \frac{1}{2}( |\uparrow_1\rangle \langle \uparrow_1|+|\downarrow_1\rangle \langle\downarrow_1| ) $$

$$S_{\rm VN_1} = -{\rm Tr}(\rho_1 \ln \rho_1) = \ln 2$$

Suppose that

$$ | \psi \rangle  =  \frac{1}{\sqrt{2}} |{\rm singlet}_{12} \rangle \otimes |\uparrow_3\rangle + \frac{1}{\sqrt{2}}|{\rm singlet}_{13} \rangle \otimes |\uparrow_2 \rangle$$

$$|\psi\rangle = \frac{1}{2} |\uparrow_1 \downarrow_2 \uparrow_3\rangle - \frac{1}{2}|\downarrow_1 \uparrow_2 \uparrow_3 \rangle + \frac{1}{2}|\uparrow_1 \uparrow_2  \downarrow_3\rangle - \frac{1}{2}|\downarrow_1 \uparrow_2  \uparrow_3 \rangle $$



$$ | \psi \rangle \langle \psi | = \frac{1}{4} |\uparrow_1 \downarrow_2 \uparrow_3 \rangle \langle \uparrow_1 \downarrow_2 \uparrow_3| - \frac{1}{4}|\downarrow_1 \uparrow_2 \uparrow_3 \rangle \langle \uparrow_1 \downarrow_2 \uparrow_3| + \frac{1}{4}|\uparrow_1 \uparrow_2  \downarrow_3 \rangle \langle \uparrow_1 \downarrow_2 \uparrow_3| - \frac{1}{4}|\downarrow_1 \uparrow_2  \uparrow_3 \rangle \langle \uparrow_1 \downarrow_2 \uparrow_3| $$

$$ -\frac{1}{4} |\uparrow_1 \downarrow_2 \uparrow_3 \rangle \langle \downarrow_1 \uparrow_2 \uparrow_3| + \frac{1}{4}|\downarrow_1 \uparrow_2 \uparrow_3 \rangle \langle \downarrow_1 \uparrow_2 \uparrow_3| - \frac{1}{4}|\uparrow_1 \uparrow_2  \downarrow_3 \rangle \langle \downarrow_1 \uparrow_2 \uparrow_3| + \frac{1}{4}|\downarrow_1 \uparrow_2  \uparrow_3 \rangle \langle \downarrow_1 \uparrow_2 \uparrow_3|$$

$$ +\frac{1}{4} |\uparrow_1 \downarrow_2 \uparrow_3 \rangle \langle \uparrow_1 \uparrow_2  \downarrow_3| - \frac{1}{4}|\downarrow_1 \uparrow_2 \uparrow_3 \rangle \langle \uparrow_1 \uparrow_2  \downarrow_3| + \frac{1}{4}|\uparrow_1 \uparrow_2  \downarrow_3 \rangle \langle \uparrow_1 \uparrow_2  \downarrow_3| - \frac{1}{4}|\downarrow_1 \uparrow_2  \uparrow_3 \rangle \langle \uparrow_1 \uparrow_2  \downarrow_3|$$

$$ -\frac{1}{4} |\uparrow_1 \downarrow_2 \uparrow_3 \rangle \langle \downarrow_1 \uparrow_2  \uparrow_3 | + \frac{1}{4}|\downarrow_1 \uparrow_2 \uparrow_3 \rangle \langle \downarrow_1 \uparrow_2  \uparrow_3 | - \frac{1}{4}|\uparrow_1 \uparrow_2  \downarrow_3 \rangle \langle \downarrow_1 \uparrow_2  \uparrow_3 | + \frac{1}{4}|\downarrow_1 \uparrow_2  \uparrow_3 \rangle \langle \downarrow_1 \uparrow_2  \uparrow_3 |$$






