from consts import *


# this module contains functions that are used to compute the ode with OBC/PBC



def x(n):
    return 2 * n * omegaF


def y(n):
    return (2 * n + 1) * omegaF


xmValsAll = [x(m) for m in range(0, N)]
ymValsAll = [y(m) for m in range(0, N)]


# mapping phi1A2
def phi1A2(deltaTau, uQ, xi):
    '''

    :param deltaTau: time interval
    :param uQ: u value at time (q+1/2)dt
    :param xi: input vector
    :return:
    '''
    for m in range(0, N):
        xi[2 * m] *= np.exp(-1j * deltaTau * (uQ + xmValsAll[m]))


# mapping phi1B2
def phi1B2(deltaTau, uQ, xi):
    '''

    :param deltaTau: time interval
    :param uQ: u value at time (q+1/2)dt
    :param xi: input vector
    :return:
    '''
    for m in range(0, N):
        xi[2 * m + 1] *= np.exp(-1j * deltaTau * (ymValsAll[m] - uQ))
    # return xi


# mapping phi1C1
def phi1C1(deltaTau, vQ, xi):
    '''

    :param deltaTau: time interval
    :param vQ: v value at time (q+1/2)dt
    :param xi: input vector
    :return:
    '''
    for m in range(0, N):
        xi[2 * m] += -1j * deltaTau * vQ * xi[2 * m + 1]
    # return xi


# mapping phi1D1
def phi1D1(deltaTau, vQ, xi):
    '''

    :param deltaTau: time interval
    :param vQ: v value at (q+1/2)dt
    :param xi: input vector
    :return:
    '''
    for m in range(0, N):
        xi[2 * m + 1] += -1j * deltaTau * vQ * xi[2 * m]


# mapping phi21
def phi21(deltaTau, wQ, xi):
    '''

    :param deltaTau: time interval
    :param wQ:  w value at time (q+1/2)dt
    :param xi: input vector
    :return:
    '''
    # j=1,2,...,N-1
    for j in range(1, N):
        xi[2 * j] += -1j * deltaTau * wQ * xi[2 * j - 1]
    # pbc
    # xi[2 * 0] += -1j * deltaTau * wQ * xi[2 * 0 - 1 + 2 * N]


# mapping phi22
def phi22(deltaTau, wQ, xi):
    '''

    :param deltaTau: time interval
    :param wQ: w value at time (q+1/2)dt
    :param xi: input  vector
    :return:
    '''
    # j=0,1,...,N-2
    for j in range(0, N - 1):
        xi[2 * j + 1] += -1j * deltaTau * wQ * xi[2 * j + 2]
    # pbc
    # xi[2 * (N - 1) + 1] += -1j * deltaTau * wQ * xi[2 * 0]


# composite mapping zeta1
def zeta1(uQ, vQ, xi):
    phi1B2(dt / 2, uQ, xi)
    phi1A2(dt / 2, uQ, xi)
    phi1C1(dt / 2, vQ, xi)
    phi1D1(dt, vQ, xi)
    phi1C1(dt / 2, vQ, xi)
    phi1B2(dt / 2, uQ, xi)
    phi1A2(dt / 2, uQ, xi)


# composite mapping zeta2
def zeta2(wQ, xi):
    phi21(dt / 4, wQ, xi)
    phi22(dt / 2, wQ, xi)
    phi21(dt / 4, wQ, xi)


# linear step in S2
def expmidtF(q, psi1):
    t = (q + 1 / 2) * dt
    uq = u(t)
    vq = v(t)
    wq = w(t)
    zeta2(wq, psi1)
    zeta1(uq, vq, psi1)
    zeta2(wq, psi1)

def S2(psiQ, q,g):
    '''

    :param psiQ: input wavefunction at time step q
    :param q: current time step
    :param g: nonlinearity
    :return: wavefunction at next time step, q+1
    '''
    # step1, nonlinear

    psiQ1 = [elem * np.exp(-1j * g * np.abs(elem) ** 2 * dt / 2) for elem in psiQ]

    # step 2, linear
    expmidtF(q, psiQ1)

    # now psiQ2 is psiQ1, because passing arguments by reference
    # step 3, nonlinear
    psiQp1 = [elem * np.exp(-1j * g * np.abs(elem) ** 2 * dt / 2) for elem in psiQ1]

    return psiQp1


def meanXAndXWd(psiQ):
    '''

    :param psiQ: wavefunction at time q
    :return: mean position at time q
    '''
    xOut = 0
    # x2Out=0
    for j in range(0, len(psiQ)):
        xOut += j * np.abs(psiQ[j]) ** 2
    #
    # for j in range(0,len(psiQ)):
    #     x2Out+=(j-xOut)**2*np.abs(psiQ[j])**2
    return xOut


def reNormalization(vec):
    '''
    normalize vector
    :param vec: input vector
    :return:
    '''
    tmp2 = 0
    for elem in vec:
        tmp2 += np.abs(elem) ** 2
    tmp = np.sqrt(tmp2)
    vec /= tmp
    return vec

