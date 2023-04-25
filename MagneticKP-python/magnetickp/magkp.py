import os,pickle
from itertools import permutations, chain
from sympy import sqrt, I, Integer,MatrixSymbol
from sympy import symbols, Mul
from sympy.matrices import Matrix
from sympy.matrices import  eye,zeros
from sympy.physics.quantum import TensorProduct

__author__ = 'Zhang Zeying'

home = os.path.dirname(os.path.abspath(__file__))+'/fd.pkl'
with open(home, 'rb') as f:
    fd = pickle.load(f)

class kpHam(object):
    r'''
    Class to construct the kp model.

    :param inputdict: The necessary input to construct the kp model, including the representation 
        matrices and rotation matrix in k space for each symmetry operators. The format of inputdict is
        {'Unitary':
        {:math:`Q_1:(D(Q_1) , R(Q_1))`,
        :math:`Q_2:(D(Q_2) , R(Q_2))`,
        ...
        },
        'Anitunitary':
        {:math:`Q_3:(D(Q_3) , R(Q_3))`,
        ...
        }}
        where :math:`Q_i` can be string type, :math:`D(Q_i)` and :math:`R(Q_i)` are the representation 
        matrix and rotation matrix in k space of :math:`Q_i`

    .. note:: 
        for the rotation matrix of anitunitary operator :math:`Q_3`, just add a minus 
        sign to the rotation matrix of the unitary part.
    
    Example usage::

        from sympy.matrices import  eye
        from sympy.matrices import Matrix
        import magnetickp
        a=magnetickp.kpHam({'Unitary':
                      {'I':(Matrix([[1,  0],[0, -1]]),Matrix([[-1,0,0],[0,-1,0],[0,0,-1]]))},
                      'Anitunitary':
                      {'IT':(Matrix([[1,  0],[0, 1]]),eye(3))}
                      })
    '''
    def __init__(self, inputdict):
        '''
        load the input
        '''
        self.input=inputdict
#        print(inputdict)
        self.debug=False
        if 'Unitary' in inputdict:
            self.dim=list(inputdict['Unitary'].values())[0][0].shape
            self.dim,_=self.dim
        else:
            self.dim=list(inputdict['Anitunitary'].values())[0][0].shape
            self.dim,_=self.dim

#        print(self.dim)
        self.__pauli_matrices=[Matrix([[1, 0],[0, 1]]),
                        Matrix([[0, 1],[1, 0]]),
                        Matrix([[0, -I],[I,  0]]),
                        Matrix([[1,  0],[0, -1]])]
        
        self.__GM = [Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
              Matrix([[0, -I, 0], [I, 0, 0], [0, 0, 0]]),
              Matrix([[0, 0, -I], [0, 0, 0], [I, 0, 0]]),
              Matrix([[0, 0, 0], [0, 0, -I], [0, I, 0]]),
              Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]),
              Matrix([[1, 0, 0], [0, -1, 0], [0, 0, 0]]),
              Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]),
              Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]),
              (1/sqrt(3)) * Matrix([[1, 0, 0], [0, 1, 0], [0, 0, -2]])]
    def __bases(self,n):
        if n==1:
            return [eye(1)]
        elif n==2:
            return self.__pauli_matrices
        elif n==3:
            return self.__GM
        elif n==4:
            return [TensorProduct(self.__pauli_matrices[i], self.__pauli_matrices[j]) for i in range(4) for j in range(4)]
        elif n==6:
            return [TensorProduct(self.__pauli_matrices[i], self.__GM[j]) for i in range(4) for j in range(9)]
        elif n==8:
            return [TensorProduct(self.__pauli_matrices[i], self.__pauli_matrices[j], self.__pauli_matrices[k]) for i in range(4) for j in range(4) for k in range(4) ]
        else:
            raise ValueError('''The Python version of MagneticKP currently supports dimension of representation matrix euqal to 1, 2, 3, 4, 6, and 8. \nFor other dimension, please use the Mathematica version of MagneticKP temporarily.''')
    def __FrobeniusInnerProduct(self,a,b):
        return ((a.H)*b).trace().nsimplify()
    def __Gf(self,g,x):
        return g*x*(g.inv())
    def __Gfc(self,g,x):
        return g*(x.C)*(g.inv())
    def getklist(self,order):
        '''
        Gererate the k bases.
        '''
        def integer_partitions(n, k, l):
            if n == 0:
                yield [0] * k
            elif k == 0:
                return
            elif n < 0:
                return
            else:
                for i in range(l, 0, -1):
                    for p in integer_partitions(n - i, k - 1, i):
                        yield [i] + p

        def getexp(expr):
            powers = {kx:0,ky:0,kz:0}
            for term in Mul.make_args(expr):
                if term.is_Pow:
                    powers[term.base] = term.exp
                else:
                    powers[term] = 1
            return(powers)
        klist = []
        kx, ky ,kz= symbols("kx ky kz")
        for partition in integer_partitions(order, 3, order):
            for perm in set(permutations(partition)):
                klist.append(kx ** perm[0] * ky ** perm[1] * kz ** perm[2])
        klist.sort(key=lambda x: [(getexp(x)[kz]),(getexp(x)[ky]),(getexp(x)[kx])])
        return(klist)



    def ISAAlg(self,Sn):
        '''
        Use ISA to obtain the common nullspace of Sn. This function is core of MagneticKP package.
        '''
        U=(Sn[0].nullspace())
        U=Matrix.hstack(*U)
        if len(Sn)==1:
            return U
        for S in Sn[1:]:
            K=Matrix.hstack(*((S*U).nullspace()))
            if K.shape==(0,0):
                return Matrix([[]])
            U=U*K
        return U
    def to_Ham(self,order,kbases,hbases,U):
        '''
        Convert the common null space of Sn to symbolic Hamiltonian
        '''
        m,n=U.shape
   #     print(m,n)
        C=[symbols(f'C_{j}_{order}') for j in range(n)]
       # print("C",C[0])
        bases=[]
        for i in kbases:
            for j in hbases:
                bases.append(i*j)
        if self.debug:
            print('base',bases)
            print('base',len(bases))
            print('base',kbases)
        ham=zeros(self.dim)
        for i in range(n):
            col=(U.col(i))
            #print(C[i])
            tem=zeros(self.dim)
            for j in range(m):
               # print(j)
                tem+=col[j]*bases[j]
        #    print(tem)

            ham+=tem*(C[i])
        if self.debug:
            print(ham)
        return ham



    def IterativelySimplify(self,korder):
        '''
        Construct the kp model for specific k-order.
        '''
        X=self.__bases(self.dim)
        basedim=Integer(korder+1)*(korder+2)/2
        if self.debug:
            print(X)
        Sn=[]
        for key in self.input.keys():
            if self.debug:
                print('123',key)
            Gg=zeros(self.dim**2)
            def getF(key):
                if fd.get(key) is not None:
                    return fd[key]
                else:
                    ValueError('''The Python version of MagneticKP currently supports rotation matrix in the SpaceGroupIrep package and k-order less than 8, i.e. the output of getRotMatOfK \nFor other cases, please use the Mathematica version of MagneticKP temporarily.''')

            if key=='Unitary':
            #   print(self.input[key].keys())
                for op in self.input[key].keys():
                    ni=-1
                    F=getF(tuple(self.input[key][op][1])+(korder,))
             #       print('F',F)
             #       print(op)
                    for j in X:
                        ni+=1
                        nj=-1
                        for i in X:
                            nj+=1
                            Gg[ni,nj]=(self.__FrobeniusInnerProduct(self.__Gf(self.input[key][op][0],i),j)/self.__FrobeniusInnerProduct(j,j))
                    S=eye((self.dim**2) *basedim)-TensorProduct(F,Gg)
                    Sn.append(S)
            elif key=='Anitunitary':
                for op in self.input[key].keys():
                    ni=-1
                    F=getF(tuple(self.input[key][op][1])+(korder,))
             #       print(self.input[key].keys())
                    for j in X:
                        nj=-1
                        ni+=1
                        for i in X:
                            nj+=1
                            Gg[ni,nj]=(self.__FrobeniusInnerProduct(self.__Gfc(self.input[key][op][0],i),j)/self.__FrobeniusInnerProduct(j,j))
                    S=eye((self.dim**2) *basedim)-TensorProduct(F,Gg)
                    Sn.append(S)
            else:
                raise ValueError('''The key of input should be either 'Unitary' or 'Anitunitary'.''')
#            print(TensorProduct(F,Gg))


#        print(Sn)
        U=self.ISAAlg(Sn)
        kbases=self.getklist(korder)
        if self.debug:
            print('kbases',korder, kbases)
        ham=self.to_Ham(korder,kbases,X,U)
        return (ham, U.shape[1])
    def getkpHam(self,order):
        '''
        Construct the kp model for specified k-order.

        :param order: can be an integer or a list of intgers. When order is
            an integer MagneticKP will calculate the kp model with k-order less 
            than or equal to order. When order is a list of intgers, MagneticKP 
            will calculate the kp model with specified k-order in the list. the
            output of getkpHam is a dictionary with format:
            {'ham': output kp Hamiltonian,
            'NumberOfParameters' : number of parameters in kp Hamiltonian,
            'korder' : considered k-order in this calculation, 
            'dim' : dimension of representation matrix in this calculation
            }

        Example usage::

            from sympy.matrices import  eye
            from sympy.matrices import Matrix
            import magnetickp
            a=magnetickp.kpHam({'Unitary':
                          {'I':(Matrix([[1,  0],[0, -1]]),Matrix([[-1,0,0],[0,-1,0],[0,0,-1]]))},
                          'Anitunitary':
                          {'IT':(Matrix([[1,  0],[0, 1]]),eye(3))}
                          })
            print(a.getkpHam(2))
            print(a.getkpHam([0,1,3]))
        '''
        ham=zeros(self.dim)
        npara=0
        if isinstance(order, int):
            korder=range(order+1)
            for i in korder:
                h,n=self.IterativelySimplify(i)
                ham+=h
                npara+=n
            return {'ham':ham,'NumberOfParameters':npara,'korder':list(korder),'dim':self.dim}
        else:
            for i in order:
                h,n=self.IterativelySimplify(i)
                ham+=h
                npara+=n
            return {'ham':ham,'NumberOfParameters':npara,'korder':list(order),'dim':self.dim}


if __name__=='__main__':
    a=kpHam({'Unitary':{'I':(Matrix([[1,  0],[0, -1]]),Matrix([[0,-1,0],[1,0,0],[0,0,1]]))},'Anitunitary':{'IT':(Matrix([[1,  0],[0, 1]]),eye(3))}})
    a.IterativelySimplify(2)

