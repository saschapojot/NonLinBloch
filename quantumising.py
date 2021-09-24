import sys
import numpy as np


class Conf(object):
    "A namespace for functions on configurations"

    def binconf(c, L):
        "binary representation of a conf integer"
        return np.binary_repr(c, L)

    def ket(conf, L):
        "quantum ket formatting with first site on the left"
        sconf = Conf.binconf(conf, L)
        return "|" + "".join([s for s in sconf[::-1]]) + ">"

    def readsite(conf, i):
        "reads quantum number on site i"
        return (conf & (1 << i)) >> i

    def conf2int(conf):
        return int("".join(str(c) for c in conf[::-1]), base=2)

    def int2conf(intconf, L):
        return [Conf.readsite(intconf, i) for i in range(L)]


class IsingHilbert(object):
    "Hilbert space for Ising in a transverse field model"

    def __init__(self, L=1):
        "description is a list containing the local hilbert space at site i"
        self.L = L
        self.hilbertsize = 2 ** self.L
        self.Confs = {i: i for i in range(self.hilbertsize)}

    def __repr__(self):
        "prints the first 30 confs of hilbert space"
        from itertools import islice
        num = min(30, self.hilbertsize)
        s = "Printing first " + str(num) + " configurations among " + str(self.hilbertsize) + "\n"
        s += ", ".join([Conf.ket(c, self.L) for c in list(islice(self.Confs.keys(), num))])
        return s

    def index(self, conf):
        "returns the index of a configuration"
        return self.Confs[conf]


def cpp_convert(x, p):
    if type(x) is complex or type(x) is np.complex64 or type(x) is np.complex128:
        return "(" + str(round(x.real, p)) + "," + str(round(x.imag, p)) + ")"
    else:
        return str(round(x, p))


class Wavefunction(np.ndarray):
    def __new__(cls, hilbert, dtype=float):
        obj = np.ndarray.__new__(cls, (hilbert.hilbertsize,), dtype)
        obj.hilbert = hilbert
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        self.hilbert = getattr(obj, "hilbert", None)

    def setvalue(self, conf, value):
        "Sets the coefficient in front of conf"
        self[self.hilbert.index(conf)] = value

    def getvalue(self, conf):
        "gets the coefficient in front of conf"
        return self[self.hilbert.index(conf)]

    def readable(self, precision=4, eps=1e-4):
        from itertools import islice
        num = min(30, self.hilbert.hilbertsize)
        L = self.hilbert.L
        s = ""
        for conf, i in islice(self.hilbert.Confs.items(), num):
            coef = self[i]
            if abs(coef) > eps:
                symb = " + "
                if type(coef) is np.float64 and np.sign(coef) == -1.0:
                    symb = " - "
                    coef = abs(coef)
                s += symb + cpp_convert(coef, precision) + Conf.ket(conf, L)
        if s[:3] == " + ":
            s = s[3:]
        return "|psi> = " + s

    def norm(self):
        return np.linalg.norm(self)

    def normalize(self):
        norm = self.norm()
        try:
            self /= norm
        except:
            print("can't normalize a vector of norm ", norm)

    def randomize(self, a=-1, b=1):
        self[:] = (b - a) * np.random.random_sample(self.size).reshape(self.shape) + a
        if self.dtype == np.complex128:
            self.imag = (b - a) * np.random.random_sample(self.size).reshape(self.shape) + a
        self.normalize()

    def measure(self, op, hermitic=True):
        "computes mean value of an hermitian operator"
        res = scalar(self, op * self)
        if abs(res.imag) > 1e-15:
            print("Warning, issue with hermiticity")
        return res.real


def scalar(psi1, psi2):
    "computes overlaps of psi1 and psi2 with complex conjugation on psi1"
    return np.vdot(psi1, psi2)


def zeroWf(hilbert, dtype=float):
    psi = Wavefunction(hilbert, dtype=dtype)
    psi[:] = 0.0
    return psi


def basisvectorWf(hilbert, index, dtype=float):
    psi = zeroWf(hilbert, dtype=dtype)
    psi[index] = 1
    return psi


def randomWf(hilbert, dtype=float):
    psi = Wavefunction(hilbert, dtype=dtype)
    psi.randomize()
    return psi


def copy(psi):
    res = Wavefunction(psi.hilbert, dtype=psi.dtype)
    res[:] = psi[:]
    return res


class Operator(object):
    def __mul__(self, psi):
        res = zeroWf(psi.hilbert, dtype=psi.dtype)
        for conf, i in psi.hilbert.Confs.items():
            coef, newconf = self.Apply(conf)
            j = psi.hilbert.Confs[newconf]
            res[j] += coef * psi[i]
        return res

    def __str__(self):
        return self.name


class Sz(Operator):
    def __init__(self, site):
        self.site = site
        self.name = "Sz_" + str(site)

    def Apply(self, conf):
        return Conf.readsite(conf, self.site) - 1 / 2, conf


class SzSz(Operator):
    def __init__(self, s1, s2):
        self.s1 = s1
        self.s2 = s2
        self.name = "Sz_" + str(s1) + "Sz_" + str(s2)

    def Apply(self, conf):
        si = Conf.readsite(conf, self.s1) - 1 / 2
        sj = Conf.readsite(conf, self.s2) - 1 / 2
        return si * sj, conf
class Sx(Operator):
    def __init__(self,site):
        self.site=site
        self.name="Sx_"+str(site)
    def Apply(self,conf):
        return 1/2, conf^(1<<self.site)


class Spinflip(Operator):
    def __init__(self,i,j):
        self.i=i
        self.j=j
        self.name="Sx_"+str(i)+"Sx_"+str(j)+"Sy_"+str(j)+"Sy_"+str(j)
    def Apply(self,conf):
        if Conf.readsite(conf,self.i)!=Conf.readsite(conf,self.j):
            return 0.5,conf^((1<<self.i)^(1<<self.j))
        else:
            return 0,conf

class GeneralOp(list):
    """
    A list of 2-element in the format [[coef, Operator]]
    corresponding to \sum_{i,j}coef_ij*Op_ij
    warning: assuming hermitian operators
    """
    def __init__(self,iterable=[]):
        list.__init__(self,iterable)
        self.having_matrix=False

    def __str__(self):
        if not len(self):
            return ""
        s="("+str(self[0][0])+")*"+str(self[0][1])
        for coef, op in self[1:]:
            s+="+("+str(coef)+")*"+str(op)
        return s
    def __iadd__(self, other):
        return GeneralOp(self[:]+other[:])

    def __add__(self, other):
        return GeneralOp(self[:]+other[:])

    def __mul__(self, psi):
        res=copy(psi)
        if self.having_matrix:
            res[:]=np.dot(self.matrix,psi)
        else:
            res[:]=0
            for coef, op in self:
                res+=coef*(op*psi)
        return res
    def createMatrix(self,hilbert,dtype=float):
        if self.having_matrix:
            return
        if hilbert.hilbertsize<2:
            return
        self.createFullMatrix(hilbert,dtype)
        self.having_matrix=True
    def createFullMatrix(self,hilbert,dtype=float):
        size=hilbert.hilbertsize
        self.matrix=np.zeros(shape=(size,size),dtype=dtype)
        for index in np.arange(size):
            self.matrix[:,index]=self*basisvectorWf(hilbert,index)

    def diagonalize(self,with_eigvec=True):
        "Assumes an hermitian operator by default"
        if not self.having_matrix:
            print("Error, don't have the matrix")
            return None
        else:
            if with_eigvec:
                self.En,self.Phin=np.linalg.eigh(self.matrix)
            else:
                self.En=np.linalg.eigvalsh(self.matrix)
    def Groundstate(self,hilbert,dtype=float):
        "computes the groundstate and returns it"
        if not self.having_matrix:
            self.createMatrix(hilbert,dtype)
        if not hasattr(self,"Phin"):
            self.diagonalize()
        gsIndex=np.argmin(self.En)
        gs=Wavefunction(hilbert,dtype=dtype)
        gs[:]=self.Phin[:,gsIndex]
        return gs


L,h=3,0.3

HamH=GeneralOp([[-h,Sx(i)] for i in range(L)])
HamIsing=GeneralOp([[-1,SzSz(i,(i+1)%L)] for i in range(L)])
Ham=HamIsing+HamH
hilbert=IsingHilbert(L=L)

psi0=Ham.Groundstate(hilbert)
print(psi0)