import numpy as np
import matplotlib.pyplot as plt



n = int(input())
xy = []
for i in range(n):
    xy.append(list(map(float,input().split())))

X = []
Y = []
for i in range(n):
    X.append(xy[i][0])
    Y.append(xy[i][1])


fig = plt.figure(figsize=(6.0, 6.0))
ax = fig.add_subplot(111)
ax.axhline(color='black',linewidth=0.5)
ax.plot(X,Y,"o",linewidth=1,color='black',label="data")



x,y = np.loadtxt("data_spline_interpolation.txt", comments='#', unpack=True)
ax.plot(x,y,"-",linewidth=1,color='black',label="spline")


#x,y = np.loadtxt("data_spline_interpolation1.txt", comments='#', unpack=True)
#ax.plot(x,y,"--",linewidth=1,color='red',label="spline1")


x,y = np.loadtxt("data_lagrange_interpolation.txt", comments='#', unpack=True)
ax.plot(x,y,"-",linewidth=1,color='blue',label="lagrange")



x,y = np.loadtxt("data_least_squares.txt", comments='#', unpack=True)
ax.plot(x,y,"-",linewidth=1,color='green',label="LSM")




ax.set_ylabel("$y$",fontsize = 11)
ax.set_xlabel("$x$",fontsize = 11)


ax.legend(loc="best")



plt.show()
fig.savefig("data_fig.png")



