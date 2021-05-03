from consts import *


inDir="./pump11/"

g=1.5

inFileName=inDir+"g"+str(g)+"wv.csv"

inDat=np.genfromtxt(inFileName,delimiter=",")
print(inDat[0,:])