import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression


g=2

inDir="./pump11/"
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
tAnn=inTVals[-1000]
yAnn=posLin[-1]-5
plt.figure(figsize=(20,20))
plt.plot(inTVals,regPosVals,color="black")
plt.plot(inTVals,posLin,color="red")
plt.annotate("y="+str(b0)+"+"+str(b1)+"t",xy=(2000,yAnn))
plt.title("g="+str(g))

outFigName=inDir+"linRegAndPos"+"g"+str(g)+".png"
plt.savefig(outFigName)