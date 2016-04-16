#Programmer: Chris Tralie
#Purpose: To create a GUI using Matplotlib for choosing surfaces of revolution
import numpy as np
import matplotlib.pyplot as plt
from PolyMesh import *
from sys import argv, exit

##A few example surfaces of revolution
def getBullet(N, NSteps):
    X = np.zeros((N, 2))
    X[:, 1] = np.pi - np.linspace(0, np.pi, N) #Go from pi down to zero
    X[:, 0] = np.sqrt(1 + np.cos(X[:, 1]))
    VPos, ITris = makeSurfaceOfRevolution(X, NSteps)
    saveOffFileExternal("Bullet.off", VPos, np.zeros(0), ITris)

def getHersheysKiss(N, NSteps):
    X = np.zeros((N, 2))
    X[:, 1] = np.pi - np.linspace(0, np.pi, N)
    X[:, 0] = 1.1+np.cos(X[:, 1])
    VPos, ITris = makeSurfaceOfRevolution(X, NSteps)
    saveOffFileExternal("HersheysKiss.off", VPos, np.zeros(0), ITris)


#Choose the surface of revolution by picking a bunch of points in the upper
#right quadrant of the XY plane
class SORChooser(object):
    def __init__(self):
        self.Points = []
        fig = plt.figure()
        self.ax = fig.add_subplot(111)
        self.ax.set_xlim((0, 10))
        self.ax.set_ylim((0, 10))
        cid = fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.drawCurve()
        self.displayAxes()
        plt.show()
        fig.canvas.mpl_disconnect(cid)
    
    def displayAxes(self):
        self.ax.set_title('Surface Of Revolution Chooser')
        self.ax.set_xlabel("X")
        self.ax.set_ylabel("Y")
    
    def drawCurve(self):
        self.ax.hold(True)
        X = np.array(self.Points)
        if X.size == 0:
            return
        self.ax.plot(X[:, 0], X[:, 1])
        self.ax.scatter(X[:, 0], X[:, 1], 20, 'r')        
            
    def onclick(self, event):
        if event.button == 1:
            self.Points.append([event.xdata, event.ydata])
        elif len(self.Points) > 0:
            self.Points.pop()
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        self.ax.clear()
        self.drawCurve()
        self.displayAxes()
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)
        plt.show()
    
    def getX(self):
        X = np.array(self.Points)
        return X

if __name__ == '__main__':
    if len(argv) < 3:
        print "Usage: python SORChooser.py <NSteps> <filename>"
        exit(0)
    s = SORChooser()
    (VPos, ITris) = makeSurfaceOfRevolution(s.getX(), int(argv[1]))
    saveOffFileExternal(argv[2], VPos, np.zeros(0), ITris)
