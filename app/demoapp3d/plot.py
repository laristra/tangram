# Simple script to read and plot the material polygons
import sys
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import matplotlib.pyplot as plt

# TODO - error checking on args
fname = sys.argv[1]
# the first two lines are diagnostics
data = np.loadtxt(fname, skiprows=2, usecols = range(25))
# to color the polygons
colors = ['y','b','g','c','m']
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# each entry is a material polygon
P=[[1,0,0],[0,1,0],[0,0,1]]
for matpoly in data:
    matid = int(matpoly[0])
    points = np.array([[matpoly[1],matpoly[2],matpoly[3]],
                       [matpoly[4],matpoly[5],matpoly[6]],
                       [matpoly[10],matpoly[11],matpoly[12]],
                       [matpoly[7],matpoly[8],matpoly[9]],
                       [matpoly[13],matpoly[14],matpoly[15]],
                       [matpoly[16],matpoly[17],matpoly[18]],
                       [matpoly[22],matpoly[23],matpoly[24]],
                       [matpoly[19],matpoly[20],matpoly[21]]])

    Z = np.zeros((8,3))
    for i in range(8): Z[i,:] = np.dot(points[i,:],P)
    r = [-1,1]
    X, Y = np.meshgrid(r, r)
    
    # list of sides' polygons of figure
    verts = [[Z[0],Z[1],Z[2],Z[3]],
             [Z[4],Z[5],Z[6],Z[7]], 
             [Z[0],Z[1],Z[5],Z[4]], 
             [Z[2],Z[3],Z[7],Z[6]], 
             [Z[1],Z[2],Z[6],Z[5]],
             [Z[4],Z[7],Z[3],Z[0]], 
             [Z[2],Z[3],Z[7],Z[6]]]
    
    # plot sides
    ax.add_collection3d(Poly3DCollection(verts, 
                                         facecolors=colors[matid%len(colors)], linewidths=1, edgecolors='k', alpha=1))


ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.savefig("%s.png" % fname.split("_")[-1], bbox_inches="tight")
plt.show()

