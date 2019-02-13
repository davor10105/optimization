EPSILON=10**-14

class Matrix():
    
    #Matrix konstruktor moze biti pozvan s listom ili tupleom, u slucaju liste, vraca se matrica sa zadanim rasporedom, u slucaju tuplea, vraca se matrica popunjena nulama dimenzija zadanih tupleom
    def __init__(self,initial):        
        if type(initial) is tuple:
            if len(initial)!=2:
                raise ValueError("Zadane dimenzije moraju biti 2-tuple")
            initial=[[0 for i in range(initial[1])] for j in range(initial[0])]
        if type(initial[0]) is not list:
            initial=[initial]
        self.matrix=[[initial[i][j] for j in range(len(initial[0]))] for i in range(len(initial))]

    #vraca dimenzije matrice u obliku (X,Y)
    def dimensions(self):
        return (len(self.matrix),len(self.matrix[0]))

    #radi kopiju matrice
    def copy(self):
        return Matrix(self.matrix)

    #zbrajanje matrica
    def __add__(self,other):
        if self.dimensions()!=other.dimensions():
            raise ValueError("Dimenzije matrice nisu odgovarajuce")
        temp=Matrix(self.matrix)
        for i in range(self.dimensions()[0]):
            for j in range(self.dimensions()[1]):
                temp.matrix[i][j]=self.matrix[i][j]+other.matrix[i][j]
        return temp

    #oduzimanje matrica
    def __sub__(self,other):
        if self.dimensions()!=other.dimensions():
            raise ValueError("Dimenzije matrice nisu odgovarajuce")
        temp=Matrix(self.matrix)
        for i in range(self.dimensions()[0]):
            for j in range(self.dimensions()[1]):
                temp.matrix[i][j]=self.matrix[i][j]-other.matrix[i][j]
        return temp

    #pomocna funkcija
    def addToTemp(temp,number,i,k):
        temp.matrix[i][k]+=number

    #mnozenje matrica
    def __rmul__(self,other):
        return self*other

    #mnozenje matrica 
    def __mul__(self,other):
        if not isinstance(other,Matrix):
            temp=self.copy()
            for i in range(self.dimensions()[0]):
                for j in range(self.dimensions()[1]):
                    temp.matrix[i][j]*=other
            return temp
    
        if self.dimensions()[1]!=other.dimensions()[0]:
            raise ValueError("Dimenzije matrice nisu odgovarajuce")
        
        temp=Matrix((self.dimensions()[0],other.dimensions()[1]))
        for i in range(self.dimensions()[0]):
                [[Matrix.addToTemp(temp,self.matrix[i][j]*other.matrix[j][k],i,k) for j in range(self.dimensions()[1])] for k in range(other.dimensions()[1])]

        return temp

    #vraca element zadan unutar indexiranja []
    def __getitem__(self,key):
        return self.matrix[key]

    #postavlja element zadan unutar indexiranja [] na item
    def __setitem__(self,key,item):
        self.matrix[key]=item

    #vraca transponiranu matricu
    def transpose(self):
        return Matrix([[self.matrix[i][j] for i in range(self.dimensions()[0])] for j in range(self.dimensions()[1])])

    #citanje matrice iz datoteke
    def read(file):
        with open(file) as f:
            temp=[]
            for line in f.readlines():
                numbers=line.rstrip().replace('\t',' ').split(' ')
                for i in range(len(numbers)):
                    numbers[i]=float(numbers[i])
                temp.append(numbers)
            return Matrix(temp)

    #pisanje matrice u datoteku
    def write(self,file):
        with open(file,'w') as f:
            f.write(str(self))

    #usporedivanje dvije matrice
    def __eq__(self,other):
        if self.dimensions()!=other.dimensions():
            return False
        X,Y=self.dimensions()
        for i in range(X):
            for j in range(Y):
                if self.matrix[i][j]!=other.matrix[i][j]:
                    return False
        return True
               
    def __str__(self):
        result=''
        for row in self.matrix:
            for element in row:
                result+=str(element)+' '
            result=result[:-1]
            result+="\n"
        return result[:-1]

    #LU dekompozicija, vraca jednu matricu koja sadrzi i L i U matricu
    def LUDecomposition(self):
        temp=self.copy()
        N,Y=self.dimensions()
        if N!=Y:
            raise ValueError

        for i in range(N-1):
            for j in range(i+1,N):
                if abs(temp.matrix[i][i])<EPSILON:
                    raise ZeroDivisionError("Stozerni element je manji od epsilon")
                temp.matrix[j][i]/=temp.matrix[i][i]
                for k in range(i+1,N):
                    temp.matrix[j][k]-=temp.matrix[i][k]*temp.matrix[j][i]
        return temp

    #LUP dekompozicija, vraca LU matricu i P matricu
    def LUPDecomposition(self):
        temp=self.copy()
        N,Y=self.dimensions()
        if N!=Y:
            raise ValueError

        P=Matrix((N,N))
        for i in range(N):
            P[i][i]=1

        for i in range(N-1):
            #nadi najveci element
            maximum=i
            for rowI in range(i+1,len(temp.matrix)):
                if abs(temp.matrix[rowI][i])>abs(temp.matrix[maximum][i]):
                    maximum=rowI
            #zamijeni redove
            t=temp.matrix[i]
            temp.matrix[i]=temp.matrix[maximum]
            temp.matrix[maximum]=t

            t=P.matrix[i]
            P.matrix[i]=P.matrix[maximum]
            P.matrix[maximum]=t

            for j in range(i+1,N):
                if abs(temp.matrix[i][i])<EPSILON:
                    raise ZeroDivisionError("Stozerni element je manji od epsilon")
                temp.matrix[j][i]/=temp.matrix[i][i]
                for k in range(i+1,N):
                    temp.matrix[j][k]-=temp.matrix[i][k]*temp.matrix[j][i]
        return (temp,P)
                    
    #supstitucija unaprijed, ako je zadana P matrica, cilajni vektor b ce se pomnoziti njome
    def forwardSupstitution(LU,b,P=None):
        N,Y=LU.dimensions()

        if P is not None:
            b=P*b
        
        y=[]
        y.append(b[0][0])
        for i in range(1,N):
            y.append(b[i][0])
            for j in range(i):
                y[i]-=LU[i][j]*y[j]

        return Matrix(y).transpose()

    #supstitucija unatrag pomocu LU matrice i vektora y
    def backwardSupstitution(LU,y):
        N,Y=LU.dimensions()

        x=[0 for i in range(N)]
        
        for i in range(N-1,-1,-1):
            if abs(LU[i][i])<EPSILON:
                raise ZeroDivisionError("Element je manji od epsilon")
            x[i]=y[i][0]
            for j in range(i+1,N):
                x[i]-=LU[i][j]*x[j]
            x[i]/=LU[i][i]
        
        return Matrix(x).transpose()

    #rjesava sustav pomocu LU dekompozicije i supstitucije unaprijed i unatrag te vraca vektor x
    def solveLU(A,b):
        LU=A.LUDecomposition()
        return Matrix.backwardSupstitution(LU,Matrix.forwardSupstitution(LU,b))

    #rjesava sustav pomocu LUP dekompozicije i supstitucije unaprijed i unatrag te vraca vektor x
    def solveLUP(A,b):
        LU,P=A.LUPDecomposition()
        return Matrix.backwardSupstitution(LU,Matrix.forwardSupstitution(LU,b,P))

    def inverse(self):
        N,Y=self.dimensions()

        if N!=Y:
            raise ValueError("Matrica nije kvadratna")

        inverseMatrix=Matrix((N,N))

        LU,P=self.LUPDecomposition()

        for i in range(N):
            iVector=[0 for x in range(N)]
            iVector[i]=1
            iVector=Matrix(iVector).transpose()

            inverseRow=Matrix.backwardSupstitution(LU,Matrix.forwardSupstitution(LU,iVector,P)).transpose()
            inverseMatrix.matrix[i]=inverseRow.matrix[0]

        return inverseMatrix.transpose()

