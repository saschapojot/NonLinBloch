import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import pandas as pd
import scipy.linalg as slin

# total lattice number(each lattice contains 2 points A and B)
N = 2 ** 10

# gaussian wavepacket parameters
# center
xc = (2 * N - 1) / 2
# xc=101
# width
sgm = 15


# parameters for coeffcient u,v w
D0 = 2
d0 =0.8
J = -1
#coef functions
def u(t):
    return D0 * np.cos(omega * t)


def v(t):
    return J + d0 * np.sin(omega * t)


def w(t):
    return J - d0 * np.sin(omega * t)

# parameters of linear part of Hamiltonian
omega = 0.03
# tilt strength
omegaF = 0.05
T = 2 * np.pi / omega
Q = 2 ** 18
tTot = 48* T
dt = tTot / Q

# init gaussian part of the wavefunction
wvFcnt = [np.exp(-(n - xc) ** 2 / (4 * sgm ** 2)) for n in range(0, 2 * N)]
# plt.figure()
# plt.plot(wvFcnt)
# plt.show()
# shape factor from one band
t0=0
k0=0
AVal=(v(t0)+w(t0))*np.cos(k0)
BVal=(v(t0)-w(t0))*np.sin(k0)
CVal=u(t0)
EVAL=np.sqrt(AVal**2+BVal**2+CVal**2)
# v1=[iB-A, C+E] without normalization const
v10=1j*BVal-AVal
v11=CVal+EVAL


for j in range(0, N):
    wvFcnt[2 * j] *= v10
    wvFcnt[2 * j + 1] *= v11
# normalization const
c02Tmp = 0
for elem in wvFcnt:
    c02Tmp += np.abs(elem) ** 2
C0 = np.sqrt(c02Tmp)

psi0 = []
for n in range(0, 2 * N):
    psi0.append(wvFcnt[n] * np.exp(1j * n * k0) / C0)

L = len(psi0)
