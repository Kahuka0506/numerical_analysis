import numpy as np
import matplotlib.pyplot as plt


fig = plt.figure(figsize=(7.0, 5.0))
ax = fig.add_subplot(111)
ax.axhline(color='black',linewidth=0.5)


x,y = np.loadtxt("./out/data_ode.txt", comments='#', unpack=True)
ax.plot(x,y,"-",linewidth=1,color='green',label="runge kutta")




ax.set_ylabel("$x$",fontsize = 11)
ax.set_xlabel("$t$",fontsize = 11)


ax.legend(loc="best")



plt.show()
fig.savefig("./out/data_ode.png")



