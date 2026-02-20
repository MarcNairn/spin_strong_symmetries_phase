#%%
import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

from qutip import *
from qutip.piqs import *

# %%# Parameters
N = 6
ntls = N
nds = num_dicke_states(ntls)
[jx, jy, jz] = jspin(N)
jp = jspin(N, "+")
jm = jp.dag()
w0 = 1
h = w0 * jz

phi = np.pi/4 #phase
#photonic parameters
nphot = 5
wc = 1
kappa = 1
ratio_g = 2
g = ratio_g/np.sqrt(N)
a = destroy(nphot)
# %%
system = Dicke(N = N)
system.hamiltonian = h 
liouv = system.liouvillian() 
system

# %%
#photonic liouvilian
h_phot = wc * a.dag() * a
c_ops_phot = [np.sqrt(kappa) * a]
liouv_phot = liouvillian(h_phot, c_ops_phot)
# %%

#identity operators
id_tls = to_super(qeye(nds))
id_phot = to_super(qeye(nphot))

# light-matter superoperator 
h_int = g * tensor(a + a.dag(), jx)
liouv_int = -1j* spre(h_int) + 1j* spost(h_int)

# total liouvillian
liouv_sum = super_tensor(liouv_phot, id_tls) + super_tensor(id_phot, liouv)
liouv_tot = liouv_sum + liouv_int
# %%
# #total operators
jz_tot = tensor(qeye(nphot), jz)
jp_tot = tensor(qeye(nphot), jp)
jm_tot = tensor(qeye(nphot), jm)
jpjm_tot = tensor(qeye(nphot), jp*jm)
nphot_tot = tensor(a.dag()*a, qeye(nds))
adag_tot = tensor(a.dag(), qeye(nds))
a_tot = tensor(a, qeye(nds))
# %%

