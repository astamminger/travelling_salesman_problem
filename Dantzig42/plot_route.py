#!/usr/bin/python
#Plot the logistic map 

from pylab import *

#load data into array
cities = loadtxt('best_route.dat')

#make axis labels
matplotlib.rc('text', usetex=True)

fig = figure()
ax1 = fig.add_subplot(111)

#plot data
print('Setting Axes Limits to: ')
length = len(cities[:,1]) - 1

maxx = max(cities[:,1])
minx = min(cities[:,1])
maxy = max(cities[:,2])
miny = min(cities[:,2])
if (maxx > maxy):
    addup = 0.01 * maxx
elif (maxy > maxx):
    addup = 0.01 * maxy
maxx = maxx + addup
minx = minx - addup
maxy = maxy + addup
miny = miny - addup

print('x_max: ',maxx)
print('x_min: ',minx)
print('y_max: ',maxy)
print('y_min: ',miny)

ax1.plot(cities[:,1], cities[:,2],c='b')
ax1.plot([cities[0,1], cities[length,1]], [cities[0,2], cities[length,2]],c='b')
ax1.scatter(cities[:,1], cities[:,2],c='g')

#set axes labels, range and title
ax1.set_xlabel(r'$x$', fontsize=18)
ax1.set_ylabel(r'$y$', fontsize=18)

ax1.set_xlim(minx,maxx)
ax1.set_ylim(miny,maxy)

#save plot and show
savefig('opt_route.svg')
show()
