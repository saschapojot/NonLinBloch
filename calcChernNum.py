from sympy import *
import numpy as np
import plotly.express as px
import matplotlib.pyplot as plt
from scipy.integrate import dblquad

J=0.0001

d0=1
D0=1

omega=0.1
T=2*np.pi/omega

def theta(p,t):
    return -2*J*d0*D0*omega*(1-np.cos(omega*t)**2*np.cos(p)**2)/\
           (4*J**2*np.cos(p)**2+4*d0**2*np.sin(omega*t)**2*np.sin(p)**2+D0**2*np.cos(omega*t)**2)**(3/2)



options={"limit":100}
rst,err=dblquad(lambda t,p:theta(p,t),0,2*np.pi,lambda p:0,lambda p:T,limit=100)
print(rst/(2*np.pi))
print(err)