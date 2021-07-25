import numpy as np
import matplotlib.pyplot as plt


fig = plt.figure(figsize=(7.0, 5.0))
ax = fig.add_subplot(111)
ax.axhline(color='black',linewidth=0.5)



x,y = np.loadtxt("./out/data_parabolic_explicit.txt", comments='#', unpack=True)
ax.plot(x,y,"-",linewidth=1,color='black',label="explicit")

x,y = np.loadtxt("./out/data_parabolic_implicit.txt", comments='#', unpack=True)
ax.plot(x,y,"-",linewidth=1,color='blue',label="implicit")


x,y = np.loadtxt("./out/data_parabolic_crank_nocolson.txt", comments='#', unpack=True)
ax.plot(x,y,"-",linewidth=1,color='red',label="crank_nocolson")




ax.set_ylabel("$u$",fontsize = 11)
ax.set_xlabel("$x$",fontsize = 11)


ax.legend(loc="best")



plt.show()
fig.savefig("./out/data_parabolic.png")



