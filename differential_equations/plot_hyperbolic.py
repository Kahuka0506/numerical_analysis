import numpy as np
import matplotlib.pyplot as plt


fig = plt.figure(figsize=(7.0, 5.0))
ax = fig.add_subplot(111)
ax.axhline(color='black',linewidth=0.5)


x,y = np.loadtxt("./out/data_hyperbollic_0.txt", comments='#', unpack=True)
ax.plot(x,y,"-",linewidth=1,color='black',label="t = 0.0")

x,y = np.loadtxt("./out/data_hyperbollic_1000.txt", comments='#', unpack=True)
ax.plot(x,y,"-",linewidth=1,color='blue',label="t = 0.2")

x,y = np.loadtxt("./out/data_hyperbollic_2000.txt", comments='#', unpack=True)
ax.plot(x,y,"-",linewidth=1,color='green',label="t = 0.4")

x,y = np.loadtxt("./out/data_hyperbollic_3000.txt", comments='#', unpack=True)
ax.plot(x,y,"-",linewidth=1,color='orange',label="t = 0.6")




ax.set_ylabel("$x$",fontsize = 11)
ax.set_xlabel("$t$",fontsize = 11)


ax.legend(loc="best")



plt.show()
fig.savefig("./out/data_hyperbollic.png")



