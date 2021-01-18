import numpy as np
import matplotlib.pyplot as plt
from FIG2_SingleLoop import func1

""" For stable solution"""
omegaArray=np.array([2,5,10,20,32])
#"""
# For unstable solution
#omegaArray=np.array([2,10,20,40,100])
omega=omegaArray[0]


# fig, ax=plt.subplots(5,1)#,sharex='col')
# ax[0].grid()
# drange=func1(omega)[0]
# u=func1(omega)[1]
# uExact=func1(omega)[2]
# l1=ax[0].plot(drange, u,drange,uExact)
# #ax[0].title.set_text(r'$\beta=$'+str(beta))
# ax[0].title.set_text(r'$\beta=0$')
DT1=func1(omega)[3]
print('Delta T= '+str(DT1))

omega=omegaArray[1]
# ax[1].grid()
# drange=func1(omega)[0]
# u=func1(omega)[1]
# uExact=func1(omega)[2]
# l2=ax[1].plot(drange, u,drange,uExact)
# ax[1].title.set_text(r'$\beta=1/12$')
DT2=func1(omega)[3]
print('Delta T= '+str(DT2))

omega=omegaArray[2]
# ax[2].grid()
# drange=func1(omega)[0]
# u=func1(omega)[1]
# uExact=func1(omega)[2]
# l3=ax[2].plot(drange, u,drange,uExact)
# ax[2].title.set_text(r'$\beta=1/6$')
DT3=func1(omega)[3]
print('Delta T= '+str(DT3))

omega=omegaArray[3]
# ax[3].grid()
# drange=func1(omega)[0]
# u=func1(omega)[1]
# uExact=func1(omega)[2]
# l4=ax[3].plot(drange, u,drange,uExact)
# ax[3].title.set_text(r'$\beta=1/4$')

DT4=func1(omega)[3]
print('Delta T= '+str(DT4))

omega=omegaArray[4]
# ax[4].grid()
# drange=func1(omega)[0]
# u=func1(omega)[1]
# uExact=func1(omega)[2]
# ax[4].plot(drange, u,drange,uExact)
# ax[4].title.set_text(r'$\beta=1/3$')
DT5=func1(omega)[3]
print('Delta T= '+str(DT5))

TArray=np.array([DT1,DT2,DT3,DT4,DT5])
plt.figure(1)
plt.plot(omegaArray,TArray,color='black', linestyle='solid', marker='o',markerfacecolor='blue', markersize=12)
plt.xlabel(r'Angular frequency ($\omega$) ')
plt.ylabel(r' Phase error ($\Delta T$)')
# ax[2].set_ylabel('Displacement (u)')
# #plt.tight_layout(pad=1, w_pad=0.5, h_pad=1.2)
# # Preventing overlaps between title
# fig.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.09)



# Labels to use in the legend for each line
"""
line_labels = ["Line A", "Line B", "Line C", "Line D"]
# Create the legend
fig.legend([l1, l2, l3, l4],     # The line objects
           labels=line_labels,   # The labels for each line
           loc="upper right",   # Position of legend
           borderaxespad=0.1,    # Small spacing around legend box
           title="Legend Title"  # Title for the legend
           )
# """
# line_labels = ["Newmark's Solution", "Exact Solution"]
# # Create the legend
# fig.legend([l1, l2],     # The line objects
#            labels=line_labels,   # The labels for each line
#            loc="center right",   # Position of legend
#            #borderaxespad=0.1,    # Small spacing around legend box
#            #title="Legend Title"  # Title for the legend
#            )

# # Adjust the scaling factor to fit your legend text completely outside the plot
# # (smaller value results in more space being made for the legend)
# fig.subplots_adjust(right=0.68)

plt.savefig('FIG_2_stable.jpg',dpi=200)