from spltFuncs import *

import inspect

outDir = "./pump1/"
gVals=[0]
tStart = datetime.now()
for g in gVals:
    dataAll=[]
    dataAll.append(psi0)
    for q in range(0,Q):
        psiCurr = dataAll[q]
        psiNext = S2(psiCurr, q,g)
        psiNext = reNormalization(psiNext)
        dataAll.append(psiNext)
    outCsvName = outDir + "g"+str(g)+"wv.csv"
    np.savetxt(outCsvName,dataAll,delimiter=",")
    xPos = []
    # wdAll=[]
    for q in range(0, Q):
        vecTmp = dataAll[q]

        xTmp = meanXAndXWd(vecTmp)
        xPos.append(xTmp)
        # wdAll.append(wdTmp)
        drift = [elem - xc for elem in xPos]
        posMax = np.max(drift)
        posMin = np.min(drift)

        # posDiff = 0.1
        # tickNum = int((posMax - posMin) / posDiff)
        # yTicks = [posMin + j * posDiff for j in range(0, tickNum + 2)]
        tAll = [dt * q for q in range(0, Q)]
        plt.figure(figsize=(20, 20))
        # plt.yticks(yTicks)
        plt.plot(tAll, drift, color="black")
        plt.xlabel("time")
        plt.ylabel("avg position")
        plt.title("g = " + str(g))
        plt.savefig(outDir + "g" + str(g) + "position.png")
        plt.close()
        # write params info

    outTxt = outDir +"g"+str(g)+ "info.txt"

    fptr = open(outTxt, "w+")
    fptr.write("Total drift = " + str(drift[-1] - drift[0]) + "\n")
    fptr.write("g=" + str(g) + "\n")
    fptr.write("omega=" + str(omega) + "\n")
    fptr.write("omegaF=" + str(omegaF) + "\n")
    fptr.write(inspect.getsource(x))
    fptr.write(inspect.getsource(y))
    fptr.write(inspect.getsource(u))
    fptr.write(inspect.getsource(v))
    fptr.write(inspect.getsource(w))
    fptr.close()
tEnd = datetime.now()
print("computation time: ", tEnd - tStart)