from spltFuncs import *


inDir="./pump11/"

g=1.5

inFileName=inDir+"g"+str(g)+"wv.csv"

inDat=np.loadtxt(inFileName,delimiter=",",dtype=complex)
rowN,colN=inDat.shape
print(rowN)
print(colN)
tAll=[dt*q for q in range(0,rowN)]
posAll=[]
for rCurr in range(0,rowN):
    posTmp=meanXAndXWd(inDat[rCurr,:])
    posAll.append(posTmp)

data={"t":tAll,"pos":posAll}

outDat=pd.DataFrame(data)
outFileName=inDir+"g"+str(g)+"position.csv"

outDat.to_csv(outFileName)
