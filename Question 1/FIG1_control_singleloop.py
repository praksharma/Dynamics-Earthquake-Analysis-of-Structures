import numpy as np
import matplotlib.pyplot as plt
from FIG1_SingleLoop import func1
betaArray=np.array([0.01,1/12,1/6,1/4,1/3])  
beta=betaArray[0]



fig, ax=plt.subplots(5,1,sharex='col')
ax[0].grid()
drange=func1(beta)[0]
u=func1(beta)[1]
uExact=func1(beta)[2]
l1=ax[0].plot(drange, u,drange,uExact)
#ax[0].title.set_text(r'$\beta=$'+str(beta))
ax[0].title.set_text(r'$\beta=0$')

print('Delta T= '+str(func1(beta)[3]))

beta=betaArray[1]
ax[1].grid()
u=func1(beta)[1]
l2=ax[1].plot(drange, u,drange,uExact)
ax[1].title.set_text(r'$\beta=1/12$')

print('Delta T= '+str(func1(beta)[3]))

beta=betaArray[2]
ax[2].grid()
u=func1(beta)[1]
l3=ax[2].plot(drange, u,drange,uExact)
ax[2].title.set_text(r'$\beta=1/6$')

print('Delta T= '+str(func1(beta)[3]))

beta=betaArray[3]
ax[3].grid()
u=func1(beta)[1]
l4=ax[3].plot(drange, u,drange,uExact)
ax[3].title.set_text(r'$\beta=1/4$')

print('Delta T= '+str(func1(beta)[3]))

beta=betaArray[4]
ax[4].grid()
u=func1(beta)[1]
ax[4].plot(drange, u,drange,uExact)
ax[4].title.set_text(r'$\beta=1/3$')


# x and y common labels https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.text.html
#fig.text(0.5, 0.02, 'Time ', ha='center')
#fig.text(0.04, 0.5, 'Displacement (u)', va='center', rotation='vertical')

#plt.subplots_adjust(left=0.5, bottom=None, right=None, top=None, wspace=None, hspace=1)

print('Delta T= '+str(func1(beta)[3]))
plt.xlabel('Time (t)')  # common label
ax[2].set_ylabel('Displacement (u)') #plt.ylabel would give label to last plot

#plt.tight_layout(pad=1, w_pad=0.5, h_pad=1.2)
# Preventing overlaps between title
fig.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.09)



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
"""
line_labels = ["Newmark's Solution", "Exact Solution"]
# Create the legend
fig.legend([l1, l2],     # The line objects
           labels=line_labels,   # The labels for each line
           loc="center right",   # Position of legend
           #borderaxespad=0.1,    # Small spacing around legend box
           #title="Legend Title"  # Title for the legend
           )

# Adjust the scaling factor to fit your legend text completely outside the plot
# (smaller value results in more space being made for the legend)
fig.subplots_adjust(right=0.68)

plt.savefig('FIG_1.jpg',dpi=200)