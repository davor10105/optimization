from Optimization import *
from utils import *

class GradientFunction():
    def __init__(self,dfList):
        self.df=dfList
        self.call_counter=0
    def __call__(self,i):
        self.call_counter+=1
        return [fun(i) for fun in self.df]

#lista ogranicenja u tupleovima zapisana donja i gornja granica
class ExplicitBoundary():
    def __init__(self,boundaries):
        self.boundaries=boundaries
    def accepted(self,x):
        newX=[]
        for i in range(len(self.boundaries)):
            if x[i]<self.boundaries[i][0]:
                newX.append(self.boundaries[i][0])
            elif x[i]>self.boundaries[i][1]:
                newX.append(self.boundaries[i][1])
            else:
                newX.append(x[i])
        return Point(tuple(newX))

#lista tupleova, (f,bool) bool true znaci da je ogranicenje jednakosti
class ImplicitBoundary():
    def __init__(self,boundaries):
        self.boundaries=boundaries
    def accepted(self,x):
        for i in range(len(self.boundaries)):
            if self.boundaries[i][1]==True:
                y=self.boundaries[i][0](x)
                if y!=0:
                    return False
            else:
                y=self.boundaries[i][0](x)
                if y<0:
                    return False
        return True

    def acceptedNonEqual(self,x):
        for i in range(len(self.boundaries)):
            if self.boundaries[i][1]==False:
                y=self.boundaries[i][0](x)
                if y<0:
                    return False
        return True

def gradientDescent(f,df,start,e=None,useGoldenRatio=True):
    if e==None:
        e=Point(tuple([0.000001 for a in start]))
    lastPoint=start
    currentPoint=start
    bestValue=f(currentPoint)
    bestCounter=0
    while(True):
        direction=df(currentPoint)
        direction=Point(tuple(direction))
        direction*=-1
        if useGoldenRatio==True:
            a,b=golden_ratio(f,currentPoint,direction)
            currentPoint=(a+b)*(1/2.)
        else:
            currentPoint+=direction

        currentValue=f(currentPoint)
        if currentValue>=bestValue:
            bestCounter+=1
        else:
            bestCounter=0
            bestValue=currentValue
        if bestCounter>=100:
            return "Postupak divergira"
        
        if checkPrecision(currentPoint,lastPoint,e):
            break
        lastPoint=currentPoint
    return currentPoint

def newtonRalphson(f,df,ddf,start,e=None,useGoldenRatio=True):
    if e==None:
        e=Point(tuple([0.000001 for a in start]))
    lastPoint=start
    currentPoint=start
    bestValue=f(currentPoint)
    bestCounter=0
    while(True):
        listH=ddf(currentPoint)
        H=[]
        n=int(len(listH)**1/2.)
        for i in range(n):
            H.append([])
            for j in range(n):
                H[i].append(listH[i*n+j])
        H=Matrix(H)

        dF=Matrix(df(currentPoint))
        dF=dF.transpose()*-1

        directionMatrix=Matrix.solveLUP(H,dF)
        direction=[]
        for i in range(directionMatrix.dimensions()[0]):
            direction.append(directionMatrix[i][0])
        direction=Point(direction)
        if useGoldenRatio==True:
            a,b=golden_ratio(f,currentPoint,direction)
            currentPoint=(a+b)*(1/2.)
        else:
            currentPoint+=direction

        currentValue=f(currentPoint)
        if currentValue>=bestValue:
            bestCounter+=1
        else:
            bestCounter=0
            bestValue=currentValue
        if bestCounter>=100:
            return "Postupak divergira"
        
        if checkPrecision(currentPoint,lastPoint,e):
            break
        lastPoint=currentPoint
    return currentPoint

def calculateCentroid(points):
    xc=Point(tuple([0 for i in range(len(points[0].point))]))
    for point in points:
        xc+=point
    return xc*(1./len(points))

import random
def box(f,start,eb,ib,e=None,alpha=1.3):
    if e==None:
        e=Point(tuple([0.000001 for a in start]))
        
    if start!=eb.accepted(start) or ib.accepted(start)==False:
        print(eb.accepted(start))
        raise ValueError("Tocka nije unutar ogranicenja")
    centroid=[]
    centroid.append(start)
    n=len(start.point)
    xc=calculateCentroid(centroid)
    for t in range(2*n-1):
        r=random.random()
        xt=[]
        try:
            for b in eb.boundaries:
                xt.append(b[0]+r*(b[1]-b[0]))
        except:
            pass
        try:
            while(ib.accepted==False):
                xt=(xt+xc)*(1/2.)
        except:
            pass
        centroid.append(Point(tuple(xt)))
        xc=calculateCentroid(centroid)

    lastCentroid=xc
    currentCentroid=xc

    while(True):
        centroid=sorted(centroid,key=lambda x:f(x),reverse=True)
        xc=calculateCentroid(centroid[1:])
        xr=(1+alpha)*xc-alpha*centroid[0]

        #eksplicitna?
        try:
            xr=eb.accepted(xr)
        except:
            pass
        #implicitna
        try:
            while(ib.accepted(xr)==False):
                xr=(xr+xc)*(1/2.)
        except:
            pass

        if f(centroid[1])<f(xr):
            xr=(xr+xc)*(1/2.)
        centroid[0]=xr

        currentCentroid=calculateCentroid(centroid)
        #print(currentCentroid)
        if checkPrecision(currentCentroid,lastCentroid,e):
            break
        lastCentroid=currentCentroid
    return currentCentroid

def ln(x):
    import math
    try:
        y=math.log(x)
    except:
        y=-1000000000000000000
    return y
def mixedTransformation(f,start,ib,e=None,t=1):
    if e==None:
        e=Point(tuple([0.000001 for a in start]))

    if ib.acceptedNonEqual(start)==False:
        newU=Function(lambda x: sum([0 if fun[1]==True else 0 if fun[0](x)>=0 else -fun[0](x) for fun in ib.boundaries]))
        start=hooke_jeeves(newU,start)
        print('Nova pocetna tocka je %s'%(start))

    newF=Function(lambda x: f(x)+sum([-1./t*ln(fun[0](x)) if fun[1]==False else t*(fun[0](x))**2  for fun in ib.boundaries]))

    currentPoint=start
    lastPoint=start
    while(True):
        currentPoint=hooke_jeeves(newF,currentPoint)

        if checkPrecision(currentPoint,lastPoint,e):
            break
        lastPoint=currentPoint
        t*=10
    return currentPoint

        
    

#xc=box(f2,Point((2,4)),ExplicitBoundary([(2,5),(2,5)]),ImplicitBoundary([(Function(lambda x:x[0]+x[1]),False)]))

#print(mixedTransformation(f2,Point((1,2)),ImplicitBoundary([(Function(lambda x:x[0]-x[1]),False)])))

#df2=GradientFunction([Function(lambda x: 2*x[0]-8),Function(lambda x:8*x[1]-16)])
#ddf2=GradientFunction([Function(lambda x: 2),Function(lambda x:0),Function(lambda x: 0),Function(lambda x:8)])
x1=Point((-1.9,2))
f1=Function(lambda x: 100*(x[1]-x[0]**2)**2+(1-x[0])**2)
df1=GradientFunction([Function(lambda x: -400*x[0]*(x[1]-x[0]**2)-2*(1-x[0])),Function(lambda x: 200*(x[1]-x[0]**2))])
ddf1=GradientFunction([Function(lambda x:-400*(x[1]-3*x[0]**2)+2),Function(lambda x:-400*x[0]),Function(lambda x:-400*x[0]),Function(lambda x:200)])

f2=Function(lambda x:(x[0]-4)**2+4*(x[1]-2)**2)
x2=Point((0.1,0.3))
df2=GradientFunction([Function(lambda x: 2*x[0]-8),Function(lambda x:8*x[1]-16)])
ddf2=GradientFunction([Function(lambda x: 2),Function(lambda x:0),Function(lambda x: 0),Function(lambda x:8)])

f3=Function(lambda x:(x[0]-2)**2+(x[1]+3)**2)
x3=Point((0,0))
df3=GradientFunction([Function(lambda x: 2*x[0]-4),Function(lambda x:2*x[1]+6)])
ddf3=GradientFunction([Function(lambda x: 2),Function(lambda x:0),Function(lambda x: 0),Function(lambda x:2)])

x4=Point((0,0))
f4=Function(lambda x:(x[0]-3)**2+(x[1])**2)
df4=GradientFunction([Function(lambda x: 2*x[0]-6),Function(lambda x:2*x[1])])
ddf4=GradientFunction([Function(lambda x: 2),Function(lambda x:0),Function(lambda x: 0),Function(lambda x:2)])
#print(newtonRalphson(f2,df2,ddf2,x2))
#print(ddf2.call_counter)
