from utils import *
import matplotlib.pyplot as plt
import numpy as np

class RungeKutta():
    def solve(A,B,x,T,tMax,printEveryNumIter=1):
        ts=[]
        xs = [[] for xi in x]

        for i in range(int(tMax / T)):
            m1=A*x+B
            m2=A*(x+T/2.*m1)+B
            m3=A*(x+T/2.*m2)+B
            m4=A*(x+T*m3)+B

            x+=T/6.*(m1+2*m2+2*m3+m4)

            if (i+1)%printEveryNumIter==0:
                PrintUtil.printVariables(i, T, x)

            ts.append((i+1)*T)
            for i in range(x.dimensions()[0]):
                xs[i].append(x[i][0])
        return (ts,xs)
class Trapez():
    def solve(A,B,x,T,tMax,printEveryNumIter=1):

        ts = []
        xs = [[] for xi in x]

        unary=Matrix(A.dimensions())
        for i in range(unary.dimensions()[0]):
            unary[i][i]=1

        for i in range(int(tMax/T)):
            R=(unary-A*T*(1./2)).inverse()*(unary+A*T*(1./2))
            S=(unary-A*T*(1./2)).inverse()*T*(1./2)*B

            x=R*x+S*2


            if (i+1)%printEveryNumIter==0:
                PrintUtil.printVariables(i, T, x)

            ts.append((i + 1) * T)
            for i in range(x.dimensions()[0]):
                xs[i].append(x[i][0])

        return (ts,xs)

class PrintUtil():
    def printVariables(i,T,x):
        varString = ""
        for k in range(x.dimensions()[0]):
            varString += "x" + str(k + 1) + ": " + str(x[k][0]) + " "
        print("t: %.2f %s" % ((i + 1) * T, varString))
    def drawVariables(ts,xs):
        i=0
        for x in xs:
            i+=1
            plt.plot(ts,x,'ro')
            plt.title("Varijabla x"+str(i))
            plt.show()

#ts,xs=RungeKutta.solve(Matrix([[0,1],[-1,0]]),Matrix((2,1)),Matrix([[0],[2]]),0.01,5)
#PrintUtil.drawVariables(ts,xs)

#ts,xs=RungeKutta.solve(Matrix([[0,1],[-200,-102]]),Matrix((2,1)),Matrix([[1],[-2]]),0.01,5)
#PrintUtil.drawVariables(ts,xs)