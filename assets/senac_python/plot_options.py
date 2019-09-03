import numpy as np
from mpl_toolkits.mplot3d import Axes3D, proj3d  
import matplotlib.pyplot as plt
import matplotlib.colors
from input_senac import RAXIS, psi0

def coolfig(xaxis,yaxis,xtext,ytext,output):
    plt.figure()
    yvectemp = np.vectorize(yaxis)
    yvec = yvectemp(xaxis)
    plt.rc('font', family='serif')
    plt.plot(xaxis, yvec)
    plt.xlabel(xtext)
    plt.ylabel(ytext)
    plt.savefig('data/%s.png' % (output))
    #plt.show()

def coolfig3D(func,xtext,ytext,ztext,output):
    xaxis,yaxis,zaxis= func
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    plt.rc('font', family='serif')
    ax.plot(xaxis,yaxis,zaxis)
    ax.set_xlim(-RAXIS[0]-RAXIS[1]-psi0, RAXIS[0]+RAXIS[1]+psi0); ax.set_ylim(-RAXIS[0]-RAXIS[1]-psi0, RAXIS[0]+RAXIS[1]+psi0); ax.set_zlim(-RAXIS[0]-RAXIS[1]-psi0, RAXIS[0]+RAXIS[1]+psi0)
    ax.set_xlabel(xtext)
    ax.set_ylabel(ytext)
    ax.set_zlabel(ztext)
    plt.savefig('data/%s.png' % (output))
    #plt.show()

def coolfigSurf3D(func,xtext,ytext,ztext,output):
    xaxis, yaxis, zaxis = func
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    plt.rc('font', family='serif')
    ax.plot_surface(xaxis,yaxis,zaxis, rstride=1, cstride=1, cmap='viridis', edgecolor='none', linewidths=0.)
    ax.set_xlim(-RAXIS[0]-RAXIS[1]-psi0, RAXIS[0]+RAXIS[1]+psi0); ax.set_ylim(-RAXIS[0]-RAXIS[1]-psi0, RAXIS[0]+RAXIS[1]+psi0); ax.set_zlim(-RAXIS[0]-RAXIS[1]-psi0, RAXIS[0]+RAXIS[1]+psi0)
    ax.set_xlabel(xtext)
    ax.set_ylabel(ytext)
    ax.set_zlabel(ztext)
    plt.savefig('data/%s.png' % (output))
    plt.savefig('data/%s.pdf' % (output), bbox_inches='tight')
    plt.show()