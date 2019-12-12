from CalVar import *
import numpy as np

def f(y, dy, x):
    if y == 0:
        return 9999999999
    else:
        return sqrt(1+dy**2)/(sqrt(-y))

def f2(y, dy, x):
    return dy * (1 + x**2 * dy)

def f3(y, dy, x):
    return sqrt(1+dy**2)

def test1():
    t = np.linspace(0,np.pi,100)
    xc = [th - np.sin(th) for th in t]
    yc = [np.cos(th) - 1 for th in t]
    a = [xc[0],yc[0]]
    b = [xc[-1],yc[-1]] 
    y, x = solve(a, b, f, 0.5, N=50)
    plt.plot(x, y, label="numerical")
    plt.plot(xc, yc,label="cycloid")
    print("for numerical F = {:.5f}".format(F(y, x, f)))
    print("for cycloid F = {:.5f}".format(F(yc, xc, f)))
    plt.legend()
    plt.show()

def test2():
    a = [1,1]
    b = [5,0.2] 
    y, x = solve(a, b, f2, 0.02, N=21, epoch=10, show=False) 
    yc = [1/xi for xi in x]
    plt.plot(x, y, label="numerical")
    plt.plot(x, yc,label="y=1/x")
    print("for numerical F = {:.5f}".format(F(y, x, f2)))
    print("for 1/x F = {:.5f}".format(F(yc, x, f2)))
    print("MSE = ", MSE(yc, y))
    plt.legend()
    plt.show()

def test3():
    a = [1,1]
    b = [2,5]
    y, x = solve(a, b, f3, 0.05,10, False)
    yc = y0(x, a, b)
    plt.plot(x, y, label="numerical")
    plt.plot(x, yc,label="y=ax+b")
    print("for numerical F = {:.2f}".format(F(y, x, f3)))
    print("for ax+b F = {:.2f}".format(F(yc, x, f3)))
    print("MSE = ", MSE(yc, y))
    plt.legend()
    plt.show()


test2()