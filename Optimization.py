class Point():
    def __init__(self,point):
        self.point=point
    def __getitem__(self,index):
        return self.point[index]
    def __add__(self,other):
        return Point(tuple([a+b for a,b in zip(self.point,other.point)]))
    def __sub__(self,other):
        return Point(tuple([a-b for a,b in zip(self.point,other.point)]))
    def __mul__(self,other):
        return Point(tuple([a*other for a in self.point]))
    def __rmul__(self,other):
        return self*other
    def norm(self):
        norm=0
        for a in self.point:
            norm+=a**2
        return norm**(1/2.)
    def __str__(self):
        return str(self.point)

class Function():
    def __init__(self,f):
        self.f=f
        self.call_counter=0
    def __call__(self,a):
        self.call_counter+=1
        return self.f(a)

def unimodal(start,h,f):
    left=start-h
    right=start+h
    center=start

    fl=f(left)
    fr=f(right)
    fc=f(center)

    if fl>fc and fr>fc:
        return (left,right)
    
    step=1
    if fl<fc:
        while (fl<fc):
            right=center
            center=left
            step*=2
            left=start-h*step

            fl=f(left)
            fc=f(center)
    else:
        while (fr<fc):
            left=center
            center=right
            step*=2
            right=start+h*step
            
            fr=f(right)
            fc=f(center)
    return (left,right)

def k():
    return 0.5*(5**(1/2.)-1)

def golden_ratio(f,start,precision):
    if type(start) is not tuple:
        interval=unimodal(start,precision,f)
    a=interval[0]
    b=interval[1]

    print(a,b)

    c=b-k()*(b-a)
    d=a+k()*(b-a)
    fc=f(c)
    fd=f(d)

    print('a=%s(%.6f) c=%s(%.6f) d=%s(%.6f) b=%s(%.6f)'%(a,f(a),c,fc,d,fd,b,f(b)))
    f.call_counter-=2

    while((b-a).norm()>precision.norm()):
        if fc<fd:
            b=d
            d=c
            c=b-k()*(b-a)

            fd=fc
            fc=f(c)
        else:
            a=c
            c=d
            d=a+k()*(b-a)

            fc=fd
            fd=f(d)
        print('a=%s(%.6f) c=%s(%.6f) d=%s(%.6f) b=%s(%.6f)'%(a,f(a),c,fc,d,fd,b,f(b)))
        f.call_counter-=2
        
    return (a,b)

def checkPrecision(currentPoint,lastPoint,precision):
    if (lastPoint-currentPoint).norm()>precision.norm():
        return False
    return True
def coordinate_search(f,start,precision=None):
    if precision==None:
        precision=Point(tuple([0.000001 for a in start]))
    lastPoint=start
    currentPoint=start
    while(True):
        for i in range(len(precision.point)):
            direction=[0 for a in start]
            direction[i]=precision[i]
            direction=Point(tuple(direction))

            a,b=golden_ratio(f,currentPoint,direction)
            currentPoint=(a+b)*(1/2.)
            
        
        if checkPrecision(currentPoint,lastPoint,precision):
            break
        lastPoint=currentPoint
    return currentPoint

def simplex_nelder_mead(f,start,precision=0.000001,step=1,alpha=1,beta=0.5,gamma=2,sigma=0.5):
    simplex_points=[start]
    for i in range(len(start.point)):
        direction=[0 for a in start]
        direction[i]=step
        direction=Point(tuple(direction))
        simplex_points.append(start+direction)
    while(True):
        minimum=f(simplex_points[0])
        min_i=0
        maximum=f(simplex_points[0])
        max_i=0
        for i in range(1,len(simplex_points)):
            value=f(simplex_points[i])
            if value<minimum:
                minimum=value
                min_i=i
            if value>maximum:
                maximum=value
                max_i=i
        centroid=Point(tuple([0 for a in start]))
        for i in range(len(simplex_points)):
            if i!=max_i:
                centroid+=simplex_points[i]
        centroid*=1./(len(simplex_points)-1)

        
        reflection=(1+alpha)*centroid-alpha*simplex_points[max_i]
        fr=f(reflection)
        if fr<minimum:
            expansion=(1-gamma)*centroid+gamma*reflection
            if f(expansion)<minimum:
                simplex_points[max_i]=expansion
            else:
                simplex_points[max_i]=reflection
        else:
            reflectionIsLarger=True
            for i in range(len(simplex_points)):
                if i==max_i:
                    continue
                if f(simplex_points[i])>fr:
                    reflectionIsLarger=False
                    break
            if reflectionIsLarger:
                if fr<maximum:
                    simplex_points[max_i]=reflection
                contraction=(1-beta)*centroid+beta*simplex_points[max_i]
                if f(contraction)<f(simplex_points[max_i]):
                    simplex_points[max_i]=contraction
                else:
                    for i in range(len(simplex_points)):
                        simplex_points[i]=sigma*(simplex_points[i]+simplex_points[min_i])
            else:
                simplex_points[max_i]=reflection
        #uvjet zaustavljanja
        fc=f(centroid)

        print('Centroid= %s(%0.6f)'%(centroid,fc))
        
        suma=0
        for point in simplex_points:
            suma+=(f(point)+fc)**2
        suma=(suma/len(simplex_points))**(1/2.)
        if suma<=precision:
            break
    retVal=Point(tuple([0 for a in start]))
    for point in simplex_points:
        retVal+=point
    return retVal*(1./len(simplex_points))

def explore(f,xp,dx):
    x=xp
    for i in range(len(xp.point)):
        fx=f(x)
        direction=[0 for j in range(len(xp.point))]
        direction[i]=dx[i]
        direction=Point(tuple(direction))
        x+=direction
        fn=f(x)
        if fn>fx:
            x-=2*direction
            fn=f(x)
            if fn>fx:
                x+=direction
    return x
        
def hooke_jeeves(f,start,dx=None,precision=None):
    if dx==None:
        dx=Point(tuple([0.5 for a in start]))
    if precision==None:
        precision=Point(tuple([0.000001 for a in start]))
    base=start
    searchstart=start
    while(True):
        new=explore(f,searchstart,dx)

        print('xb=%s(%.6f) xp=%s(%.6f) xn=%s(%.6f)'%(base,f(base),searchstart,f(searchstart),new,f(new)))
        f.call_counter-=3
        
        if f(new)<f(base):
            searchstart=2*new-base
            base=new
        else:
            dx*=1/2.
            searchstart=base
        if dx.norm()<=precision.norm():
            return base
        
f1=Function(lambda x: 100*(x[1]-x[0]**2)**2+(1-x[0])**2)
x1=Point((1,1))
f2=Function(lambda x: (x[0]-4)**2+4*(x[1]-2)**2)
x2=Point((0.1,0.3))
f3=Function(lambda x:  sum([(x[i]-(i+1))**2 for i in range(len(x.point))]))
x3=Point((0,0,0,0,0))
f4=Function(lambda x: abs((x[0]-x[1])*(x[0]+x[1]))+(x[0]**2+x[1]**2)**(1/2.))
x4=Point((5.1,1.1))
import math
f5=Function(lambda x: 0.5+(math.sin(sum([i**2 for i in x.point])**(1/2.))**2-0.5)/(1+0.001*sum([i**2 for i in x.point])))
x5=Point((0,0,0))
#print(coordinate_search(f4,x4))
#print(f4.call_counter)
