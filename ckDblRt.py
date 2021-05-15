from consts import *
import plotly.graph_objects as go
#This script checks the double roots of eigenequation





def gm(t):
    return v(t) + w(t)


def theta(t):
    return v(t) - w(t)


kGridNum = 200
tGridNum = 200
deltaK = 2 * np.pi / kGridNum
deltaT = T / tGridNum

tValsAll = [deltaT * q for q in range(0, tGridNum + 1)]
kValsAll = [deltaK * j for j in range(0, kGridNum + 1)]

gVals = [1]

for g in gVals:
    tKEAll = []
    for tTmp in tValsAll:
        for kTmp in kValsAll:
            pTmp = omegaF * tTmp - kTmp
            uTmp = u(tTmp)
            gammaTmp = gm(tTmp)
            thetaTmp = theta(tTmp)
            # nonlinear
            coefTmp = [1,
                       -3 * g,
                       13 / 4 * g ** 2 - uTmp ** 2 - gammaTmp ** 2 * np.cos(pTmp) ** 2 - thetaTmp ** 2 * np.sin(
                           pTmp) ** 2,
                       -3 / 2 * g ** 3 + uTmp ** 2 * g + 2 * g * gammaTmp ** 2 * np.cos(
                           pTmp) ** 2 + 2 * g * thetaTmp ** 2 * np.sin(pTmp) ** 2,
                       1 / 4 * g ** 4 - 1 / 4 * uTmp ** 2 * g ** 2 - g ** 2 * gammaTmp ** 2 * np.cos(
                           pTmp) ** 2 - g ** 2 * thetaTmp ** 2 * np.sin(pTmp) ** 2]

            rtsTmp = np.roots(coefTmp)
            if sum(np.isreal(rtsTmp))>2:
                for rtTmp in rtsTmp:
                    if np.isreal(rtTmp):
                        tKEAll.append([tTmp, kTmp, np.real(rtTmp)])

    # finished solving eqn
    # np.savetxt("rootsg"+str(g)+".txt",tKEAll)
    tToPlot = []
    kToPlot = []
    eToplot = []
    for elem in tKEAll:
        tToPlot.append(elem[0])
        kToPlot.append(elem[1])
        eToplot.append(elem[2])

    fig = plt.figure(figsize=(20, 20))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(tToPlot, kToPlot, eToplot, marker=".", c="black")
    ax.set_xlabel("t")
    ax.set_ylabel("k")
    ax.set_zlabel("E")
    plt.title("g=" + str(g))

    outDir = "./pbcSpec2/"
    outDirFile = outDir +"dblRoots"+ "g" + str(g)
    plt.savefig(outDirFile + ".png")

    plt.close()
    solData = {"t": tToPlot, "k": kToPlot, "energy": eToplot}
    outSolCsvName = outDir +"dblRoots" +"g" + str(g) + ".csv"
    outDt = pd.DataFrame(solData)
    outDt.to_csv(outSolCsvName)

    figGo = go.Figure(data=[go.Scatter3d(
            x=tToPlot,
            y=kToPlot,
            z=eToplot,
            mode='markers',
            marker=dict(
                size=1,
                color="black",  # set color to an array/list of desired values
                # choose a colorscale
                opacity=0.8
            ),

        )])
    figGo.update_layout(scene=dict(xaxis_title="t",
                                       yaxis_title="k",
                                       zaxis_title="E"),
                            title_text="g=" + str(g))
    figGo.write_html(outDirFile +"dblRoots"+ ".html")

