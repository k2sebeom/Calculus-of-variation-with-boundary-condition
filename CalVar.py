import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, sin, pi, cos

#Helper function for numerical differentiation
def diff(y, x):
    dy = []
    for i in range(1, len(x)):
        dy.append(((y[i]- y[i-1])/(x[i]-x[i-1])))
    dy.insert(0,y[0])
    return np.array(dy)

#Helper functino for numerical definite integration
def integrate(y, x):
    I = 0
    for i in range(1, len(x)):
        I += y[i] * (x[i]-x[i-1])
    return I

#functional F based on the generator f
def F(y, x, f):
    dy = diff(y, x)
    fg = []
    for i in range(len(x)):
      fg.append(f(y[i], dy[i], x[i]))
    return integrate(fg, x)

#generates a line between point a and point b
def y0(x, a, b):
    return (a[1]-b[1])/(a[0]-b[0])*(x-a[0])+a[1]

#n th order term of the fourier series
def fourier(n, xs):
    return np.array([sin(n*pi/(xs[-1]-xs[0])*(x - xs[0])) for x in xs])

#generates a finite fourier seires based on the sequence Cn
def fseries(Cn, x):
    y = np.zeros(len(x))
    for n in range(len(Cn)):
        y += Cn[n] * fourier(n+1, x)
    return y

#reterns dF/dc
def dFdc(f, n, y1, x):
    dc = 0.0001
    y2 = y1 + dc  * fourier(n, x)
    F1 = F(y1, x, f)
    F2 = F(y2, x, f)
    return (F2- F1)/dc

#reterns n dimensional dF/dCn
def dFdCn(f, y, x, Cn):
    dF = np.zeros(len(Cn))
    for n in range(len(Cn)):
        dF[n] = dFdc(f, n+1, y, x)
    return dF
    
def norm(L):
    return np.sum(np.square(L))/len(L)

def solve2(a, b, f, lr, o=5, show=False):
    x = np.linspace(a[0],b[0],100)
    L = x[-1]-x[0]
    y = y0(x, a, b)
    Cn = np.zeros(o)
    dF = np.ones(len(Cn))
    i = 0
    Lr = np.zeros(o)
    for n in range(o):
        Lr[n] = lr/(n+1)**2
    y1 = y + fseries(Cn, x)
    print("="*5+"solving"+"="*5)
    while norm(dF) > 0.000001:
        i+=1
        if show and (i%50 == 1):
            plt.plot(x, y1)
        dF = dFdCn(f, y1, x, Cn)
        Cn = Cn - dF * Lr
        y1 = y + fseries(Cn, x)
        print("|dF/dC|^2 = {:.10f}".format(norm(dF)), end='\r')
        if norm(dF) > 1000:
            print("Learning rate is too large. Choose smaller learning rate.")
            break
    solution = "solution = y0(x)"
    for i in range(len(Cn)):
        solution += " + " + "{:.2f}".format(Cn[i]) + "sin[" + "{:.2f}".format((i+1)*pi/L) + "(x - " + str(x[0]) + ")]"
    print('\n'+"="*17)
    print(solution)
    return y1, x

def optimize(n, f, a, b, x, y):
    Fs = []
    fou = fourier(n, x)
    cs = np.linspace(-1,1,100)
    for c in cs:
        yc = y + c * fou
        Fs.append(F(yc, x, f))
    plt.plot(cs, Fs)
    plt.xlabel("c")
    plt.ylabel("F[y]")
    plt.show()

def descend(y, x, n, f, lr):
    c = 0
    fou = fourier(n, x)
    y1 = y + c * fou
    dFdc = 10
    dc = 0.0001
    lr = lr/n**2
    i = 0
    while(abs(dFdc) > 0.0001):
        i+=1
        F1 = F(y1, x, f)
        y2 = y + (c+dc)*fou
        F2 = F(y2, x, f)
        dFdc = (F2-F1)/dc
        print("optimizing order {} term: loss = {:.5f}".format(n, abs(dFdc)), end='\r')
        c -= dFdc * lr
        y1 = y + c * fou
    return y1, c

def MSE(y1, y2):
    mse = 0
    for i in range(len(y1)):
        mse += (y1[i] - y2[i])**2
    return sqrt(mse)

def solve(a, b, f, lr, N=10, epoch=3, show=False):
    x = np.linspace(a[0],b[0],2000)
    y = y0(x, a, b)
    cs = []
    print("==solving==")
    for j in range(epoch):
        print("Epoch {} / {}".format(j+1, epoch))
        for i in range(1,o):
            y, c = descend(y, x, i, f, lr)
            cs.append(c)
        print('\n'+"Result: F = {:.4f}".format(j+1, epoch,F(y, x, f)))
        if show:
            plt.plot(x, y)
    print("="*10)
    return y, x
