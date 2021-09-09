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
sgm = 50


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
T = 2* np.pi / omega
Q = 2 ** 14
tTot = 3*T
dt = tTot / Q

# init gaussian part of the wavefunction
# wvFcnt = [np.exp(-(n - xc) ** 2 / (4 * sgm ** 2)) for n in range(0, 2 * N)]

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
uA=v10
uB=v11


# for j in range(0, N):
#     wvFcnt[2 * j] *= v10
#     wvFcnt[2 * j + 1] *= v11
# # normalization const
# c02Tmp = 0
# for elem in wvFcnt:
#     c02Tmp += np.abs(elem) ** 2
# C0 = np.sqrt(c02Tmp)
#
# psi0 = []
# for n in range(0, 2 * N):
#     psi0.append(wvFcnt[n] * np.exp(1j * n * k0) / C0)



#########################################
def l2Norm(vec):
    tmp=0
    for j in range(0,len(vec)):
        tmp+=np.abs(vec[j])**2
    return np.sqrt(tmp)
j0=xc
kIndAll=list(range(0,2*N))
kValsAll=[2*np.pi*kIndTmp/(2*N) for kIndTmp in kIndAll]
fValsAll=[np.exp(-(kValsAll[kIndTmp]+k0)**2*sgm**2+1j*j0*(kValsAll[kIndTmp]+k0)) for kIndTmp in kIndAll]
psi0=np.zeros(2*N,dtype=complex)

for n in range(0,N):
    #psi2n
    psi2nVal=0
    for kIndTmp in kIndAll:
        psi2nVal+=fValsAll[kIndTmp]*np.exp(-1j*kValsAll[kIndTmp]*2*n)

    psi0[2*n]=psi2nVal*uA
    #psi2n+1
    psi2np1Val=0
    for kIndTmp in kIndAll:
        psi2np1Val+=fValsAll[kIndTmp]*np.exp(-1j*kValsAll[kIndTmp]*(2*n+1))
    psi0[2*n+1]=psi2np1Val*uB

psi0/=l2Norm(psi0)
print(np.max(np.abs(psi0))**2/np.min(np.abs(psi0))**2)
plt.figure()
plt.plot(range(0,2*N),np.abs(psi0)**2,color="black")
plt.show()


