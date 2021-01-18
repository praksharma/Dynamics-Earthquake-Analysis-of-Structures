###########################################################################
# program: generalized_eigen.py
# author: Tom Irvine  
# Email: tomirvine@aol.com
# version: 1.5
# date: September 11, 2013
# description:  generalized eigenvalue function
#
###########################################################################

from scipy.linalg import eig

import numpy as np

###########################################################################

def generalized_eigen(K,M,p):
    """
    K=stiffness matrix
    M=mass matrix
    
    p=1 print natural frequencies
    p=2 print natural frequencies and mode shapes
    
    There is no print out for other values of p

    FN = natural frequencies
    MS = mode shapes

    ndof = number of natural frequencies
    
    """    

    nlength=len(np.atleast_1d(M)) 
    
    print (" nlength = %d " %nlength) 

    if(nlength==1):   
        
        ndof=1
        
        FN=np.zeros((1,1),'f')
        omega=np.sqrt(K/M)
        oma=np.array(omega)
        FN[0]=oma/(2*np.pi)
        MS=1/np.sqrt(M)    
        
        i=0        
        
        print (" ")
        print (" i   fn(Hz)")        
        print ("%d. %8.4g" %(i,FN[i]))
        
        print (" ")
        print (" Modeshape")
        print (MS)         
        
    else:
        
        (L,V) = eig(K,M)

        omega=np.sqrt(L)
        oma=np.array(omega)
        fn=oma/(2*np.pi)
        fn=abs(fn)

        order=fn.ravel().argsort()

        ndof=len(L)
        mf=ndof

        NN=np.zeros(ndof,'f')
        FN=np.zeros(ndof,'f')

        for i in range (0,ndof):
            NN[i]=float(i)
            FN[i]=fn[order[i]]
        
#    
#   Mass Normalize Eigenvectors
#
        QQQ=np.dot(V.T,np.dot(M,V))
    
        for i in range (0,mf):
            nf=np.sqrt(QQQ[i,i])
            for j in range (0,mf):
                V[j,i]/=nf      
#
#   Sort Eigenvectors
#
        MS=np.zeros((mf,mf),'f')
#
        for i in range (0,mf):
            MS[0:mf,i]=V[0:mf,order[i]]          

        mfs=ndof
        if(mfs>100):
            mfs=100
        
        if(p==1 or p==2):    
            print (" ")
            print (" i   fn(Hz)")
            for i in range (0,mfs):
                print ("%d. %8.4g" %(i,FN[i]))
            
        if(p==2):
            print (" ")
            print (" Modeshapes")
            print (MS)    
        
    return ndof,FN,MS    