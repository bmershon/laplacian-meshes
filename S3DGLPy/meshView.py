#Based off of http://wiki.wxpython.org/GLCanvas
#Lots of help from http://wiki.wxpython.org/Getting%20Started
from OpenGL.GL import *
from OpenGL.arrays import vbo
import wx
from wx import glcanvas

from Primitives3D import *
from PolyMesh import *
from LaplacianMesh import *
from Cameras3D import *
from MeshCanvas import *
from struct import *
from sys import exit, argv
import random
import numpy as np
import scipy.io as sio
from pylab import cm
import os
import math
import time
from time import sleep
from pylab import cm
import matplotlib.pyplot as plt

#GUI States
(STATE_NORMAL, STATE_SAVEROTATIONS, STATE_SAVELIGHTING, STATE_CHOOSELAPLACEVERTICES) = (0, 1, 2, 3)
#Laplacian substates
(SUBSTATE_NONE, CHOOSELAPLACE_WAITING, CHOOSELAPLACE_PICKVERTEX) = (0, 1, 2)

class MeshViewerCanvas(BasicMeshCanvas):
    def __init__(self, parent):
        super(MeshViewerCanvas, self).__init__(parent)
        self.GUIState = STATE_NORMAL
        self.GUISubstate = SUBSTATE_NONE
        
        #State variables for laplacian mesh operations
        self.laplacianConstraints = {} #Elements will be key-value pairs (idx, Point3D(new position))
        self.laplaceCurrentIdx = -1
    
    def displayMeshFacesCheckbox(self, evt):
        self.displayMeshFaces = evt.Checked()
        self.Refresh()

    def displayMeshEdgesCheckbox(self, evt):
        self.displayMeshEdges = evt.Checked()
        self.Refresh()

    def displayMeshVerticesCheckbox(self, evt):
        self.displayMeshVertices = evt.Checked()
        self.Refresh()

    def displayVertexNormalsCheckbox(self, evt):
        self.displayVertexNormals = evt.Checked()
        self.Refresh()

    def displayFaceNormalsCheckbox(self, evt):
        self.displayFaceNormals = evt.Checked()
        self.Refresh()

    def useLightingCheckbox(self, evt):
        self.useLighting = evt.Checked()
        self.needsDisplayUpdate = True
        self.Refresh()
    
    def useTextureCheckbox(self, evt):
        self.useTexture = evt.Checked()
        self.needsDisplayUpdate = True
        self.Refresh()    
    
    #######Laplacian mesh menu handles
    def doLaplacianMeshSelectVertices(self, evt):
        if self.mesh:
            self.mesh.updateIndexDisplayList()
            self.GUIState = STATE_CHOOSELAPLACEVERTICES
            self.GUISubstate = CHOOSELAPLACE_WAITING
            self.Refresh()
    
    def clearLaplacianMeshSelection(self, evt):
        self.laplacianConstraints.clear()
        self.Refresh()
    
    def doLaplacianSolveWithConstraints(self, evt):
        anchorWeights = 1e8
        anchors = np.zeros((len(self.laplacianConstraints), 3))
        i = 0
        anchorsIdx = []
        for anchor in self.laplacianConstraints:
            anchorsIdx.append(anchor)
            anchors[i, :] = self.laplacianConstraints[anchor]
            i += 1
        
        #IGL Cotangent weights
        (L, M_inv, solver, deltaCoords) = makeLaplacianMatrixSolverIGLSoft(self.mesh.VPos, self.mesh.ITris, anchorsIdx, anchorWeights)
        self.mesh.VPos = solveLaplacianMatrixIGLSoft(solver, L, M_inv, deltaCoords, anchorsIdx, anchors, anchorWeights)
        
#        #My umbrella weights
#        L = makeLaplacianMatrixUmbrellaWeights(self.mesh.VPos, self.mesh.ITris, anchorsIdx, anchorWeights)
#        deltaCoords = L.dot(self.mesh.VPos)[0:self.mesh.VPos.shape[0], :]
#        self.mesh.VPos = np.array(solveLaplacianMatrix(L, deltaCoords, anchors, anchorWeights), dtype=np.float32)
        
        sio.savemat("anchors.mat", {'deltaCoords':deltaCoords, 'anchors':anchors, 'anchorsIdx':np.array(anchorsIdx)})
        self.mesh.needsDisplayUpdate = True
        self.mesh.updateIndexDisplayList()
        self.Refresh()

    def computeMeanCurvatures(self, evt):
        (L, solver, deltaCoords) = makeLaplacianMatrixSolverIGL(self.mesh.VPos, self.mesh.ITris, anchorsIdx, anchorWeights)
        #Color vertices to be equal to mean curvature
        curvs = np.sqrt(np.sum(deltaCoords**2, 1))
        curvs = curvs - np.min(curvs)
        curvs = curvs / np.max(curvs)
        cmConvert = cm.get_cmap('jet')
        self.mesh.VColors = cmConvert(curvs)[:, 0:3]
        self.mesh.needsDisplayUpdate = True

    def drawMeshStandard(self):
        glEnable(GL_LIGHTING)
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, [0.8, 0.8, 0.8, 1.0])
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, [0.2, 0.2, 0.2, 1.0])
        glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, 64)
        self.camera.gotoCameraFrame()
        glLightfv(GL_LIGHT0, GL_POSITION, np.array([0, 0, 0, 1]))
        self.mesh.renderGL(self.displayMeshEdges, self.displayMeshVertices, self.displayMeshFaces, self.displayVertexNormals, self.displayFaceNormals, self.useLighting, self.useTexture)

    def repaint(self):
        self.setupPerspectiveMatrix()
        
        glClearColor(0.0, 0.0, 0.0, 0.0)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        
        if self.GUIState == STATE_NORMAL and self.mesh:
            self.drawMeshStandard()        
        elif self.GUIState == STATE_CHOOSELAPLACEVERTICES:
            if self.GUISubstate == CHOOSELAPLACE_WAITING:
                glDisable(GL_LIGHTING)
                self.camera.gotoCameraFrame()
                glPointSize(10)
                glBegin(GL_POINTS)
                for idx in self.laplacianConstraints:
                    P = self.mesh.VPos[self.mesh.vertices[idx].ID, :]
                    glColor3f(0, 1, 0)
                    glVertex3f(P[0], P[1], P[2])
                    P = self.laplacianConstraints[idx]
                    if idx == self.laplaceCurrentIdx:
                        glColor3f(1, 0, 0)
                    else:
                        glColor3f(0, 0, 1)
                    glVertex3f(P[0], P[1], P[2])
                glEnd()
                glColor3f(1, 1, 0)
                glBegin(GL_LINES)
                for idx in self.laplacianConstraints:
                    P1 = self.mesh.VPos[self.mesh.vertices[idx].ID, :]
                    P2 = self.laplacianConstraints[idx]
                    glVertex3f(P1[0], P1[1], P1[2])
                    glVertex3f(P2[0], P2[1], P2[2])
                glEnd()
                self.drawMeshStandard()
                
            elif self.GUISubstate == CHOOSELAPLACE_PICKVERTEX:
                glClearColor(0.0, 0.0, 0.0, 0.0)
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
                self.mesh.renderGLIndices()
                pixel = glReadPixels(self.MousePos[0], self.MousePos[1], 1, 1, GL_RGBA, GL_UNSIGNED_BYTE)
                [R, G, B, A] = [int(pixel.encode("hex")[i*2:(i+1)*2], 16) for i in range(4)]
                idx = extractFromRGBA(R, G, B, 0) - 1
                if idx >= 0 and idx < len(self.mesh.vertices):
                    if idx in self.laplacianConstraints:
                        #De-select if it's already selected
                        self.laplaceCurrentIdx = -1
                        self.laplacianConstraints.pop(idx, None)
                    else:
                        self.laplacianConstraints[idx] = np.array(self.mesh.VPos[self.mesh.vertices[idx].ID, :])
                        self.laplaceCurrentIdx = idx
                self.GUISubstate = CHOOSELAPLACE_WAITING
                self.Refresh()
        self.SwapBuffers()

    def MouseDown(self, evt):
        state = wx.GetMouseState()
        if self.GUIState == STATE_CHOOSELAPLACEVERTICES:
            if state.ShiftDown():
                #Pick vertex for laplacian mesh constraints
                self.GUISubstate = CHOOSELAPLACE_PICKVERTEX
        x, y = evt.GetPosition()
        self.CaptureMouse()
        self.handleMouseStuff(x, y)
        self.Refresh()

    def MouseMotion(self, evt):
        state = wx.GetMouseState()
        x, y = evt.GetPosition()
        [lastX, lastY] = self.MousePos
        self.handleMouseStuff(x, y)
        dX = self.MousePos[0] - lastX
        dY = self.MousePos[1] - lastY
        if evt.Dragging():
            idx = self.laplaceCurrentIdx
            if self.GUIState == STATE_CHOOSELAPLACEVERTICES and state.ControlDown() and self.laplaceCurrentIdx in self.laplacianConstraints:
                #Move up laplacian mesh constraint based on where the user drags
                #the mouse
                t = self.camera.towards
                u = self.camera.up
                r = np.cross(t, u)
                P0 = self.laplacianConstraints[idx]
                #Construct a plane going through the point which is parallel to the
                #viewing plane
                plane = Plane3D(P0, t)
                #Construct a ray through the pixel where the user is clicking
                tanFOV = math.tan(self.camera.yfov/2)
                scaleX = tanFOV*(self.MousePos[0] - self.size.x/2)/(self.size.x/2)
                scaleY = tanFOV*(self.MousePos[1] - self.size.y/2)/(self.size.y/2)
                V = t + scaleX*r + scaleY*u
                ray = Ray3D(self.camera.eye, V)
                self.laplacianConstraints[idx] = ray.intersectPlane(plane)[1]
            else:
                #Translate/rotate shape
                if evt.MiddleIsDown():
                    self.camera.translate(dX, dY)
                elif evt.RightIsDown():
                    self.camera.zoom(-dY)#Want to zoom in as the mouse goes up
                elif evt.LeftIsDown():
                    self.camera.orbitLeftRight(dX)
                    self.camera.orbitUpDown(dY)
        self.Refresh()

class MeshViewerFrame(wx.Frame):
    (ID_LOADDATASET, ID_SAVEDATASET, ID_SAVEDATASETMETERS, ID_SAVESCREENSHOT, ID_CONNECTEDCOMPONENTS, ID_SPLITFACES, ID_TRUNCATE, ID_FILLHOLES, ID_GEODESICDISTANCES, ID_PRST, ID_INTERPOLATECOLORS, ID_SAVEROTATINGSCREENSOTS, ID_SAVELIGHTINGSCREENSHOTS, ID_SELECTLAPLACEVERTICES, ID_CLEARLAPLACEVERTICES, ID_SOLVEWITHCONSTRAINTS, ID_MEMBRANEWITHCONSTRAINTS, ID_GETHKS, ID_GETHEATFLOW) = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)
    
    def __init__(self, parent, id, title, pos=DEFAULT_POS, size=DEFAULT_SIZE, style=wx.DEFAULT_FRAME_STYLE, name = 'GLWindow'):
        style = style | wx.NO_FULL_REPAINT_ON_RESIZE
        super(MeshViewerFrame, self).__init__(parent, id, title, pos, size, style, name)
        #Initialize the menu
        self.CreateStatusBar()
        
        self.size = size
        self.pos = pos
        print "MeshViewerFrameSize = %s, pos = %s"%(self.size, self.pos)
        self.glcanvas = MeshViewerCanvas(self)
        
        #####File menu
        filemenu = wx.Menu()
        menuOpenMesh = filemenu.Append(MeshViewerFrame.ID_LOADDATASET, "&Load Mesh","Load a polygon mesh")
        self.Bind(wx.EVT_MENU, self.OnLoadMesh, menuOpenMesh)
        menuSaveMesh = filemenu.Append(MeshViewerFrame.ID_SAVEDATASET, "&Save Mesh", "Save the edited polygon mesh")
        self.Bind(wx.EVT_MENU, self.OnSaveMesh, menuSaveMesh)
        menuSaveScreenshot = filemenu.Append(MeshViewerFrame.ID_SAVESCREENSHOT, "&Save Screenshot", "Save a screenshot of the GL Canvas")
        self.Bind(wx.EVT_MENU, self.OnSaveScreenshot, menuSaveScreenshot)
        
        menuExit = filemenu.Append(wx.ID_EXIT,"E&xit"," Terminate the program")
        self.Bind(wx.EVT_MENU, self.OnExit, menuExit)
        
        #####Laplacian Mesh Menu
        laplacianMenu = wx.Menu()
        menuSelectLaplaceVertices = laplacianMenu.Append(MeshViewerFrame.ID_SELECTLAPLACEVERTICES, "&Select Laplace Vertices", "Select Laplace Vertices")
        self.Bind(wx.EVT_MENU, self.glcanvas.doLaplacianMeshSelectVertices, menuSelectLaplaceVertices)
        
        menuClearLaplaceVertices = laplacianMenu.Append(MeshViewerFrame.ID_CLEARLAPLACEVERTICES, "&Clear vertex selection", "Clear Vertex Selection")
        self.Bind(wx.EVT_MENU, self.glcanvas.clearLaplacianMeshSelection, menuClearLaplaceVertices)
        
        menuSolveWithConstraints = laplacianMenu.Append(MeshViewerFrame.ID_SOLVEWITHCONSTRAINTS, "&Solve with Constraints", "Solve with Constraints")
        self.Bind(wx.EVT_MENU, self.glcanvas.doLaplacianSolveWithConstraints, menuSolveWithConstraints)
        
        # Creating the menubar.
        menuBar = wx.MenuBar()
        menuBar.Append(filemenu,"&File") # Adding the "filemenu" to the MenuBar
        menuBar.Append(laplacianMenu,"&MeshLaplacian") # Adding the "filemenu" to the MenuBar
        self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.
        
        self.rightPanel = wx.BoxSizer(wx.VERTICAL)
        
        #Buttons to go to a default view
        viewPanel = wx.BoxSizer(wx.HORIZONTAL)
        topViewButton = wx.Button(self, -1, "Top")
        self.Bind(wx.EVT_BUTTON, self.glcanvas.viewFromTop, topViewButton)
        viewPanel.Add(topViewButton, 0, wx.EXPAND)
        sideViewButton = wx.Button(self, -1, "Side")
        self.Bind(wx.EVT_BUTTON, self.glcanvas.viewFromSide, sideViewButton)
        viewPanel.Add(sideViewButton, 0, wx.EXPAND)
        frontViewButton = wx.Button(self, -1, "Front")
        self.Bind(wx.EVT_BUTTON, self.glcanvas.viewFromFront, frontViewButton)
        viewPanel.Add(frontViewButton, 0, wx.EXPAND)
        self.rightPanel.Add(wx.StaticText(self, label="Views"), 0, wx.EXPAND)
        self.rightPanel.Add(viewPanel, 0, wx.EXPAND)
        
        #Checkboxes for displaying data
        self.displayMeshFacesCheckbox = wx.CheckBox(self, label = "Display Mesh Faces")
        self.displayMeshFacesCheckbox.SetValue(True)
        self.Bind(wx.EVT_CHECKBOX, self.glcanvas.displayMeshFacesCheckbox, self.displayMeshFacesCheckbox)
        self.rightPanel.Add(self.displayMeshFacesCheckbox, 0, wx.EXPAND)
        self.displayMeshEdgesCheckbox = wx.CheckBox(self, label = "Display Mesh Edges")
        self.displayMeshEdgesCheckbox.SetValue(False)
        self.Bind(wx.EVT_CHECKBOX, self.glcanvas.displayMeshEdgesCheckbox, self.displayMeshEdgesCheckbox)
        self.rightPanel.Add(self.displayMeshEdgesCheckbox, 0, wx.EXPAND)
        self.displayMeshVerticesCheckbox = wx.CheckBox(self, label = "Display Mesh Points")
        self.displayMeshVerticesCheckbox.SetValue(True)
        self.Bind(wx.EVT_CHECKBOX, self.glcanvas.displayMeshVerticesCheckbox, self.displayMeshVerticesCheckbox)
        self.rightPanel.Add(self.displayMeshVerticesCheckbox, 0, wx.EXPAND)
        self.displayVertexNormalsCheckbox = wx.CheckBox(self, label = "Display Vertex Normals")
        self.displayVertexNormalsCheckbox.SetValue(False)
        self.Bind(wx.EVT_CHECKBOX, self.glcanvas.displayVertexNormalsCheckbox, self.displayVertexNormalsCheckbox)
        self.rightPanel.Add(self.displayVertexNormalsCheckbox, 0, wx.EXPAND)
        self.displayFaceNormalsCheckbox = wx.CheckBox(self, label = "Display Face Normals")
        self.displayFaceNormalsCheckbox.SetValue(False)
        self.Bind(wx.EVT_CHECKBOX, self.glcanvas.displayFaceNormalsCheckbox, self.displayFaceNormalsCheckbox)
        self.rightPanel.Add(self.displayFaceNormalsCheckbox, 0, wx.EXPAND)
        self.useLightingCheckbox = wx.CheckBox(self, label = "Use Lighting")
        self.useLightingCheckbox.SetValue(True)
        self.Bind(wx.EVT_CHECKBOX, self.glcanvas.useLightingCheckbox, self.useLightingCheckbox)
        self.rightPanel.Add(self.useLightingCheckbox, 0, wx.EXPAND)
        self.useTextureCheckbox = wx.CheckBox(self, label = "Use Texture")
        self.useTextureCheckbox.SetValue(False)
        self.Bind(wx.EVT_CHECKBOX, self.glcanvas.useTextureCheckbox, self.useTextureCheckbox)
        self.rightPanel.Add(self.useTextureCheckbox, 0, wx.EXPAND)

        #Finally add the two main panels to the sizer        
        self.sizer = wx.BoxSizer(wx.HORIZONTAL)
        #cubecanvas = CubeCanvas(self)
        #self.sizer.Add(cubecanvas, 2, wx.EXPAND)
        self.sizer.Add(self.glcanvas, 2, wx.EXPAND)
        self.sizer.Add(self.rightPanel, 0, wx.EXPAND)
        
        self.SetSizer(self.sizer)
        self.Layout()
        self.glcanvas.Show()
    
    def OnLoadMesh(self, evt):
        dlg = wx.FileDialog(self, "Choose a file", ".", "", "OFF files (*.off)|*.off|TOFF files (*.toff)|*.toff|OBJ files (*.obj)|*.obj", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            filepath = os.path.join(dirname, filename)
            print dirname
            self.glcanvas.mesh = PolyMesh()
            print "Loading mesh %s..."%filename
            self.glcanvas.mesh.loadFile(filepath)
            self.glcanvas.meshCentroid = self.glcanvas.mesh.getCentroid()
            self.glcanvas.meshPrincipalAxes = self.glcanvas.mesh.getPrincipalAxes()
            print "Finished loading mesh"
            print self.glcanvas.mesh
            self.glcanvas.initMeshBBox()
            self.glcanvas.Refresh()
        dlg.Destroy()
        return

    def OnSaveMesh(self, evt):
        dlg = wx.FileDialog(self, "Choose a file", ".", "", "*", wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            filepath = os.path.join(dirname, filename)
            self.glcanvas.mesh.saveFile(filepath, True)
            self.glcanvas.Refresh()
        dlg.Destroy()
        return     
        
    def OnSaveScreenshot(self, evt):
        dlg = wx.FileDialog(self, "Choose a file", ".", "", "*", wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            filepath = os.path.join(dirname, filename)
            saveImageGL(self.glcanvas, filepath)
        dlg.Destroy()
        return

    def OnExit(self, evt):
        self.Close(True)
        return

class MeshViewer(object):
    def __init__(self, filename = None, ts = False, sp = "", ra = 0):
        app = wx.App()
        frame = MeshViewerFrame(None, -1, 'MeshViewer')
        frame.Show(True)
        app.MainLoop()
        app.Destroy()

if __name__ == '__main__':
    viewer = MeshViewer()
