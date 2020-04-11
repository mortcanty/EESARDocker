# ;+
# ; Return digamma function (derivative of gamma function divided by
# ; gamma function)
# ; :Params:
# ;    x, in, required
# ;       Argument for digamma function
# ;    eps, in, optional, default=1e-12
# ;       Precision for result
# ; :History:
# ;    15 Apr 2008 Written, Anthony Smith
# ;-

import numpy as np
import math

def ajs_digamma(x,eps=1e-12):
#  Euler-Mascheroni constant
    gamma = 0.57721566490153286060651209008240243104215933593992
    psi = -gamma
    delta_psi = 1.0
    n = 1
    while abs(delta_psi)>eps:
        delta_psi = (x - 1) / (n * (n + x - 1))
        psi += delta_psi
        n += 1
    return psi

Llu = np.zeros((800,3))
digamma = np.zeros(800)
for L in range(10,800):
    digamma[L] = ajs_digamma(L/10.0) 
    if L%100 == 0:
        print( '%i '%L )
for L in range(10,800):
    Llu[L,0] = -digamma[L] + math.log(L/10.0)
for L in range(20,800):
    Llu[L,1] = -(2*digamma[L-10] + 1.0/(L/10.0-1.0)) + 2*math.log(L/10.0)    
for L in range(30,800):
    Llu[L,2] = -(3*digamma[L-20] + 2.0/(L/10.0-2.0) + 1.0/(L/10.0-1.0)) + 3*math.log(L/10.0) 
with open('/home/mort/Dokumente/lookup.txt','w') as f:
    f.write(Llu)    