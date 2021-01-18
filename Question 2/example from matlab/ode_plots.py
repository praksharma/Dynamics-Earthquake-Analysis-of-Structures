##########################################################################
# program: ode_plots.py
# author: Tom Irvine
# version: 1.1
# date: April 30, 2012
#
##########################################################################

import matplotlib.pyplot as plt

from tompy import enter_int

from numpy import zeros

##########################################################################

def ode_plots(fignum,nt,ndof,t,d,v,a,iu):
    """
    fignum = figure number
    nt = number of time points
    ndof = number of degrees-of-freedom
    t = time
    d = displacement
    v = velocity
    a = acceleration
    iu = units:  1=English 2=metric
    """

    plt.figure(fignum)
    fignum=fignum+1

    for j in range(0,ndof):
        mass='dof '+str(j+1)
        plt.plot(t, d[:,j], label=mass)

    plt.grid(True)
    plt.title('Displacement')
    plt.legend(loc="upper right")
    plt.xlabel('Time(sec)')

    if(iu==1):
        plt.ylabel('Disp(in)')
    else:
        plt.ylabel('Disp(m)')


    plt.figure(fignum)
    fignum=fignum+1

    for j in range(0,ndof):
        mass='dof '+str(j+1)
        plt.plot(t, v[:,j], label=mass)

    plt.grid(True)
    plt.title('Velocity')
    plt.legend(loc="upper right")
    plt.xlabel('Time(sec)')

    if(iu==1):
        plt.ylabel('Vel(in/sec)')
    else:
        plt.ylabel('Vel(m/sec)')


    plt.figure(fignum)
    fignum=fignum+1

    for j in range(0,ndof):
        mass='dof '+str(j+1)
        plt.plot(t, a[:,j], label=mass)

    plt.grid(True)
    plt.title('Acceleration')
    plt.legend(loc="upper right")
    plt.xlabel('Time(sec)')

    plt.ylabel('Accel(G)')


    print " "
    print " Plot relative displacement?"
    print "  1=yes  2=no "

    ird=enter_int();

    while(ird==1):

        plt.figure(fignum)
        fignum=fignum+1

        print "  "
        print " Enter first dof number "
        ia=enter_int()
        print "  "
        print " Enter second dof number "
        ib=enter_int()

        rd=zeros(nt,float)
        for j in range (0,nt):
            rd[j]=d[j,(ia-1)]-d[j,(ib-1)]

        plt.plot(t,rd)

        plt.grid(True)

        title_string= ' Relative Displacement  dof '+str(ia)+'-'+str(ib)

        plt.title(title_string)
        plt.xlabel('Time(sec)')

        if(iu==1):
            plt.ylabel('Disp(in)')
        else:
            plt.ylabel('Disp(m)')

        print " "
        print " Plot another relative displacement?"
        print "  1=yes  2=no "

        ird=enter_int()

    print " "
    print " view plots"

##########################################################################

def ode_plots_sdof(fignum,nt,ndof,t,d,v,a,iu):
    """
    fignum = figure number
    nt = number of time points
    ndof = number of degrees-of-freedom
    t = time
    d = displacement
    v = velocity
    a = acceleration
    iu = units:  1=English 2=metric
    """

    plt.figure(fignum)
    fignum=fignum+1

    plt.plot(t,d)

    plt.grid(True)
    plt.title('Displacement')
    plt.xlabel('Time(sec)')

    if(iu==1):
        plt.ylabel('Disp(in)')
    else:
        plt.ylabel('Disp(m)')


    plt.figure(fignum)
    fignum=fignum+1


    plt.plot(t, v)


    plt.grid(True)
    plt.title('Velocity')
    plt.xlabel('Time(sec)')

    if(iu==1):
        plt.ylabel('Vel(in/sec)')
    else:
        plt.ylabel('Vel(m/sec)')


    plt.figure(fignum)
    fignum=fignum+1


    plt.plot(t, a)

    plt.grid(True)
    plt.title('Acceleration')
    plt.xlabel('Time(sec)')

    plt.ylabel('Accel(G)')

    print " "
    print " view plots"
