##########################################################################
# program: Newmark.py
# author: Tom Irvine
# version: 1.4
# date: May 1, 2012
# description:  solution of a system of second-order ODEs for a dynamic
#               system via the Newmark-beta method
#
##########################################################################

from numpy import zeros,dot,atleast_1d

from scipy import linalg


################################################################################

def Newmark_coefficients(dt):
    alpha=0.25
    beta=0.5

    a0=1/(alpha*(dt**2))
    a1=beta/(alpha*dt)
    a2=1/(alpha*dt)
    a3=(1/(2*alpha))-1
    a4=(beta/alpha)-1
    a5=(dt/2)*((beta/alpha)-2)
    a6=dt*(1-beta)
    a7=beta*dt

    return a0,a1,a2,a3,a4,a5,a6,a7

################################################################################

def Newmark_initialize(ndof,a0,a1,M,C,K,NT,DI,VI):

    nlength=len(atleast_1d(M)) 

    KH=zeros((ndof,ndof),float)

    KH=K+a0*M+a1*C

    if(nlength==1):

        U=zeros(NT,float)
        Ud=zeros(NT,float)
        Udd=zeros(NT,float)

        U[0]=DI
        Ud[0]=VI
            
    else:

        U=zeros((ndof,NT),float)
        Ud=zeros((ndof,NT),float)
        Udd=zeros((ndof,NT),float)

        U[:,0]=DI
        Ud[:,0]=VI

    return U,Ud,Udd,KH

################################################################################

def Newmark_free(M,C,K,VI,DI,dt,NT,ndof):
    """
    input

    M = mass matrix
    C = damping matrix
    K = stiffness matrix

    VI = initial velocity
    DI = initial displacement

      dt = time step
      NT = number of time points
    ndof = number of degrees of freedom

    output

      U = displacement
     Ud = velocity
    Udd = acceleration

    """

    a0,a1,a2,a3,a4,a5,a6,a7=Newmark_coefficients(dt)

    U,Ud,Udd,KH=Newmark_initialize(ndof,a0,a1,M,C,K,NT,DI,VI)
    
    nlength=len(atleast_1d(M))     

    print "M"
    print M

    print "C"
    print C

    print "K"
    print K

    print "KH"
    print KH
    
    if(nlength==1):
 
        for i in range (1,NT):

            V1=(a1*U[i-1]+a4*Ud[i-1]+a5*Udd[i-1])
            V2=(a0*U[i-1]+a2*Ud[i-1]+a3*Udd[i-1])

            CV=dot(C,V1)
            MA=dot(M,V2)

            FH=MA+CV

#  solve for displacements

            Un=FH/KH

            Uddn=a0*(Un-U[i-1])-a2*Ud[i-1]-a3*Udd[i-1]
            Udn=Ud[i-1]+a6*Udd[i-1]+a7*Uddn

            U[i]=Un[:]
            Ud[i]=Udn[:]
            Udd[i]=Uddn[:]
            
    else:
        
  

        for i in range (1,NT):

            V1=(a1*U[:,i-1]+a4*Ud[:,i-1]+a5*Udd[:,i-1])
            V2=(a0*U[:,i-1]+a2*Ud[:,i-1]+a3*Udd[:,i-1])

            CV=dot(C,V1)
            MA=dot(M,V2)

            FH=MA+CV

#  solve for displacements

            if(ndof>1):
                Un = linalg.solve(KH, FH)
            else:
                Un=FH/KH

            Uddn=a0*(Un-U[:,i-1])-a2*Ud[:,i-1]-a3*Udd[:,i-1]
            Udn=Ud[:,i-1]+a6*Udd[:,i-1]+a7*Uddn

            U[:,i]=Un[:]
            Ud[:,i]=Udn[:]
            Udd[:,i]=Uddn[:]

    return U.T,Ud.T,Udd.T

################################################################################

def Newmark_force(M,C,K,FFI,force_dof,VI,DI,dt,NT,ndof):
    """
    input

    M = mass matrix
    C = damping matrix
    K = stiffness matrix

    FFI = interpolated force matrix
    force_dof = connects forces with dofs

    VI = initial velocity
    DI = initial displacement

      dt = time step
      NT = number of time points
    ndof = number of degrees of freedom

    output

      U = displacement
     Ud = velocity
    Udd = acceleration

    """

    a0,a1,a2,a3,a4,a5,a6,a7=Newmark_coefficients(dt)

    U,Ud,Udd,KH=Newmark_initialize(ndof,a0,a1,M,C,K,NT,DI,VI)


    for i in range (1,NT):

        V1=(a1*U[:,i-1]+a4*Ud[:,i-1]+a5*Udd[:,i-1])
        V2=(a0*U[:,i-1]+a2*Ud[:,i-1]+a3*Udd[:,i-1])

        CV=dot(C,V1)
        MA=dot(M,V2)


#  apply forces

        F=zeros(ndof,float)

        for j in range (0,ndof):

            j_index=force_dof[j]

            if(j_index!=-999):

                F[j]=FFI[i,j_index]


        FH=F+MA+CV

#  solve for displacements

        if(ndof>1):
            Un = linalg.solve(KH, FH)
        else:
            Un=FH/KH

        Uddn=a0*(Un-U[:,i-1])-a2*Ud[:,i-1]-a3*Udd[:,i-1]
        Udn=Ud[:,i-1]+a6*Udd[:,i-1]+a7*Uddn


        U[:,i]=Un
        Ud[:,i]=Udn
        Udd[:,i]=Uddn

    return U.T,Ud.T,Udd.T

################################################################################

def Newmark_force_modal(M,C,K,FFI,force_dof,VI,DI,dt,NT,ndof,MS):
    """
    input

    M = mass matrix
    C = damping matrix
    K = stiffness matrix

    FFI = interpolated force matrix
    force_dof = connects forces with dofs

    VI = initial velocity
    DI = initial displacement

      dt = time step
      NT = number of time points
    ndof = number of degrees of freedom

    MS = mode shapes

    output

      U = displacement
     Ud = velocity
    Udd = acceleration

    """
    
    nlength=len(atleast_1d(M))

    a0,a1,a2,a3,a4,a5,a6,a7=Newmark_coefficients(dt)


    if(nlength==1):

        U=zeros(NT,float)
        Ud=zeros(NT,float)
        Udd=zeros(NT,float)

        KH=K+a0*M+a1*C

        U[0]=DI
        Ud[0]=VI
        

        for i in range (1,NT):

#  apply force
            
            for j in range (0,ndof):

                j_index=force_dof[j]

                if(j_index!=-999):

                    F=FFI[i,j_index]

                    nF=MS.T*F

                V1=(a1*U[i-1]+a4*Ud[i-1]+a5*Udd[i-1])
                V2=(a0*U[i-1]+a2*Ud[i-1]+a3*Udd[i-1])

                CV=C*V1
                MA=M*V2

                FH=nF+MA+CV

#  solve for displacements

                Un=FH/KH

                Uddn=a0*(Un-U[i-1])-a2*Ud[i-1]-a3*Udd[i-1]
                Udn=Ud[i-1]+a6*Udd[i-1]+a7*Uddn

                U[i]=Un
                Ud[i]=Udn
                Udd[i]=Uddn
        
        
    else:

        U=zeros((ndof,NT),float)
        Ud=zeros((ndof,NT),float)
        Udd=zeros((ndof,NT),float)

        KH=zeros(ndof,float)
        FH=zeros(ndof,float)
        


        for j in range (0,ndof):
            KH[j]=K[j]+a0*M[j]+a1*C[j]

            U[j,0]=DI[j]
            Ud[j,0]=VI[j]

  

        for i in range (1,NT):

#  apply forces

            F=zeros(ndof,float)

            for j in range (0,ndof):

                j_index=force_dof[j]

                if(j_index!=-999):

                    F[j]=FFI[i,j_index]

                    nF=dot(MS.T,F)

            for j in range (0,ndof):

                V1=(a1*U[j,i-1]+a4*Ud[j,i-1]+a5*Udd[j,i-1])
                V2=(a0*U[j,i-1]+a2*Ud[j,i-1]+a3*Udd[j,i-1])

                CV=C[j]*V1
                MA=M[j]*V2

                FH[j]=nF[j]+MA+CV

#  solve for displacements

                Un=FH[j]/KH[j]

                Uddn=a0*(Un-U[j,i-1])-a2*Ud[j,i-1]-a3*Udd[j,i-1]
                Udn=Ud[j,i-1]+a6*Udd[j,i-1]+a7*Uddn

                U[j,i]=Un
                Ud[j,i]=Udn
                Udd[j,i]=Uddn


    return U.T,Ud.T,Udd.T
