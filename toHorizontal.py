import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

#this script removes the linear tendline

g=0
inDir="./pump1/"
inFileName=inDir+"g"+str(g)+"position.csv"
inCsv=pd.read_csv(inFileName)
inTVals=np.array(inCsv.iloc[:,1])
inPosVals=np.array(inCsv.iloc[:,2])
posStart=inPosVals[0]
regPosVals=inPosVals-posStart
regTVals=inTVals.reshape((-1,1))

md=LinearRegression().fit(regTVals,regPosVals)
b0=md.intercept_
b1=md.coef_[0]
posLin=md.predict(regTVals)

horizontalVals=[regPosVals[n]-posLin[n] for n in range(0,len(posLin))]

fSpaceHorVals=np.fft.fft(horizontalVals)

fSpaceAbs=np.abs(np.fft.fftshift(fSpaceHorVals))
frqs=np.fft.fftshift(np.fft.fftfreq(n=len(fSpaceAbs)))

indDecreasing=sorted(range(len(fSpaceAbs)),key=lambda k: fSpaceAbs[k])

