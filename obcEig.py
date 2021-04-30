from consts import *
from scipy.sparse import diags
import scipy.sparse.linalg as sslin
import scipy.linalg as slin
from scipy.sparse import csr_matrix

xmValsAll = [x(m) for m in range(0, N)]
ymValsAll = [y(m) for m in range(0, N)]

NN = N


def H0(t):
    uVal = u(t)
    vVal = v(t)
    wVal = w(t)

    d0 = []
    for j in range(0, NN):
        d0.append(uVal + xmValsAll[j])
        d0.append(-uVal + ymValsAll[j])
    d1 = []
    for j in range(0, NN - 1):
        d1.append(vVal)
        d1.append(wVal)
    d1.append(vVal)

    offSet = [-1, 0, 1]
    elems = np.array([d1, d0, d1])
    rst = diags(elems, offSet).toarray()
    return rst


def H1(u, g):
    '''

    :param u: state vector, length=2*NN
    :param g: nonlinearity
    :return: nonlinear part of Hamiltonian
    '''
    d0 = []
    for j in range(0, len(u)):
        d0.append(g * np.abs(u[j]) ** 2)
    rst = diags([d0], [0]).toarray()
    return rst


def R(u, g, H0Val):
    '''

    :param u: normalized vector
    :param g: nonlinearity
    :param H0Val: linear part of Hamiltonian
    :return: Rayleigh quotient
    '''

    rst0 = np.transpose(np.conj(u)).dot(H0Val).dot(u)
    rst1 = sum([g * np.abs(elem) ** 4 for elem in u])

    return rst0 + rst1


def nextPsi(H0Val, H1Val, ECurr, psiCurr):
    '''

    :param H0Val:  current H0
    :param H1Val: current H1
    :param ECurr: current E
    :param psiCurr: current vector psi
    :return: psi in next iteration
    '''

    d0 = np.ones(2 * NN) * ECurr
    diagTmp = diags([d0], [0]).toarray()
    matTmp = csr_matrix(H0Val + H1Val - diagTmp)
    yNext = sslin.spsolve(matTmp, psiCurr)
    normYNext = slin.norm(yNext, 2)
    yNext /= normYNext
    return yNext


M = 100
eps = 1e-15

def checkDist(eVec):
    '''

    :param eVec: input energy eigenvalue vector
    :return: True for break
    '''
    if(len(eVec)==1):
        return False
    else:
        dist=np.abs(eVec[-1]-eVec[-2])
        if dist<eps:
            return True
        else:
            return False


def energyAtt(t, g):
    '''

    :param t: time
    :return: energy spectrum at t
    '''
    H0MatCurr = H0(t)  # dense

    H0ValSp = csr_matrix(H0MatCurr)  # sparse
    EOut = []
    vecOut = []


    eigsAll, vecsAll = sslin.eigsh(H0ValSp,k=10,which="SM")
    rowN,colN=vecsAll.shape


    for j in range(0, colN):
        EjAllTmp = []

        vecj = vecsAll[:, j]  # init vec
        EjInit = R(vecj, g, H0MatCurr)
        EjAllTmp.append(EjInit)
        for m in range(0, M):
            tfVal=checkDist(EjAllTmp)
            if tfVal:
                break

            H1Curr = H1(vecj, g)
            vecj = nextPsi(H0MatCurr, H1Curr, EjAllTmp[m], vecj)
            ENext = R(vecj, g, H0MatCurr)
            EjAllTmp.append(ENext)

        EOut.append(EjAllTmp[-1])
        vecOut.append(vecj)

    return EOut, vecOut


def checkNorm(t,u,eVal,g):
    h0Val=H0(t)
    h1Val=H1(u,g)
    h=h0Val+h1Val
    err=np.abs(slin.norm(h.dot(u)-eVal*u,2))
    return err



startTime=datetime.now()
# for j in range(0,len(eVals)):
#     errTmp=checkNorm(t,vecVals[j],eVals[j],g)
#     if errMax<errTmp:
#         errMax=errTmp
# print(errMax)

TGrid=100
deltaT=T/TGrid
eValsAll=[]
vecAll=[]
qRange=range(0,TGrid+1)
tAll=[qVal*deltaT for qVal in qRange]
g=1.2
for tCurr in tAll:
    eCurr,vecCurr=energyAtt(tCurr,g)
    eValsAll.append(eCurr)


pltTVals=[]
pltEVals=[]
for q in qRange:
    for elem in eValsAll[q]:
        pltTVals.append(tAll[q])
        pltEVals.append(elem)

plt.figure(figsize=(20,20))
plt.plot(pltTVals,pltEVals,color="black")
plt.xlabel("time")
plt.ylabel("energy")
plt.title("g="+str(g))

outDir="./spec9/"
plt.savefig(outDir + "g" + str(g) + "position.png")

plt.close()
# write params info

outTxt = outDir +"g"+str(g)+ "info.txt"

fptr = open(outTxt, "w+")

fptr.write("g=" + str(g) + "\n")
fptr.write("omega=" + str(omega) + "\n")
fptr.write("omegaF=" + str(omegaF) + "\n")
fptr.write("d0="+str(d0)+"\n")
fptr.write("D0="+str(D0)+"\n")
fptr.write("J="+str(J)+"\n")
fptr.write("k0="+str(k0)+"\n")
fptr.write("sigma="+str(sgm)+"\n")
fptr.write(inspect.getsource(x))
fptr.write(inspect.getsource(y))
fptr.write(inspect.getsource(u))
fptr.write(inspect.getsource(v))
fptr.write(inspect.getsource(w))
fptr.close()

endTime=datetime.now()
print("time :", endTime-startTime)








