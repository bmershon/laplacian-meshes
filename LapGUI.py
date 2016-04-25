#Based off of http://wiki.wxpython.org/GLCanvas
#Lots of help from http://wiki.wxpython.org/Getting%20Started
import sys
sys.path.append("S3DGLPy")
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
from ColorTextureTools import *

#LAPLACIAN MESH CONSTANTS
SPECTRUM_K = 20
LOWPASS_K = 20
HEAT_K = 200
HKS_K = 200
HKS_T = 20

#GUI States
(STATE_NORMAL, STATE_CHOOSELAPLACEVERTICES, STATE_CHOOSECOLORVERTICES, STATE_ANIMATEHEAT) = (0, 1, 2, 3)
#Laplacian substates
(SUBSTATE_NONE, CHOOSELAPLACE_WAITING, CHOOSELAPLACE_PICKVERTEX) = (0, 1, 2)
#Color picking substates
(COLORPICK_NONE, COLORPICK_WAITING, COLORPICK_PICKVERTEX, COLORPICK_PICKCOLOR) = (0, 1, 2, 3)

class MeshViewerCanvas(BasicMeshCanvas):
    def clearAllSelections(self):
        #State variables for laplacian mesh operations
        self.laplacianConstraints = {} #Elements will be key-value pairs (idx, Point3D(new position))
        self.laplaceCurrentIdx = -1
        self.laplacianSelections = [] #Stores an ordered list of the selections (used for flattening)
        
        self.colorChoices = {} #Elements will be key-value paris (idx, np.array([R, G, B]))
        self.colorCurrentIdx = -1    

    def __init__(self, parent):
        super(MeshViewerCanvas, self).__init__(parent)
        self.GUIState = STATE_NORMAL
        self.GUISubstate = SUBSTATE_NONE
        
        self.clearAllSelections()
        
        #State variables for heat, etc
        (self.eigvalues, self.eigvectors) = (np.array([]), np.array([]))
        self.heatIdx = 0
        self.heat_ts = np.linspace(0, 1000, 100)
        
        #State variables for color picking
        self.colorPickTexID = None
        self.colorPosVBO = None
        self.colorColorVBO = None
        
        self.Bind(wx.EVT_KEY_DOWN, self.onKeyPress)
        self.Bind(wx.EVT_KEY_UP, self.onKeyRelease)
        self.pressingC = False
    
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
    
    def onKeyPress(self, evt):
        if evt.GetUnicodeKey() == 67:
            self.pressingC = True
    
    def onKeyRelease(self, evt):
        if evt.GetUnicodeKey() == 67:
            self.pressingC = False
    
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
    
    def getAnchors(self):
        anchors = np.zeros((len(self.laplacianConstraints), 3))
        i = 0
        anchorsIdx = []
        for anchor in self.laplacianConstraints:
            anchorsIdx.append(anchor)
            anchors[i, :] = self.laplacianConstraints[anchor]
            i += 1
        anchorsIdx = np.array(anchorsIdx)
        return (anchors, anchorsIdx)
    
    def getSelectedColors(self):
        colors = np.zeros((len(self.colorChoices), 3))
        colorsIdx = []
        i = 0
        for idx in self.colorChoices:
            colorsIdx.append(idx)
            colors[i, :] = self.colorChoices[idx]
            i += 1
        return (colors, colorsIdx)
    
    def doLaplacianSolveWithConstraints(self, evt):
        (anchors, anchorsIdx) = self.getAnchors()        
        solveLaplacianMesh(self.mesh, anchors, anchorsIdx)
        self.mesh.needsDisplayUpdate = True
        self.mesh.updateIndexDisplayList()
        self.Refresh()
    
    def estimateCurvature(self, evt):
        #Color vertices to be equal to mean curvature
        curvs = estimateMeanCurvature(self.mesh)
        curvs = curvs - np.min(curvs)
        curvs = curvs / np.max(curvs)
        cmConvert = cm.get_cmap('jet')
        self.mesh.VColors = cmConvert(curvs)[:, 0:3]
        self.mesh.needsDisplayUpdate = True
    
    def doLaplacianSmooth(self, evt):
        doLaplacianSmooth(self.mesh)
        self.mesh.needsDisplayUpdate = True
        self.Refresh()

    def doLaplacianSharpen(self, evt):
        doLaplacianSharpen(self.mesh)
        self.mesh.needsDisplayUpdate = True
        self.Refresh()
    
    def doMinimalSurface(self, evt):
        (anchors, anchorsIdx) = self.getAnchors()        
        makeMinimalSurface(self.mesh, anchors, anchorsIdx)
        #solveLaplacianMesh(self.mesh, anchors, anchorsIdx)
        self.mesh.needsDisplayUpdate = True
        self.mesh.updateIndexDisplayList()
        self.Refresh()

    def getSpectrum(self, evt):
        cmap = plt.get_cmap('jet')
        (VPos, ITris) = (self.mesh.VPos, self.mesh.ITris)
        K = SPECTRUM_K
        if (K >= self.mesh.VPos.shape[0]):
            K = self.mesh.VPos.shape[0]-1
        (lam, U) = getLaplacianSpectrum(self.mesh, K)
        for k in range(U.shape[1]):
            Y = np.zeros(U.shape[0])
            if k > 0:
                Y = U[:, k]
                Y = Y - np.min(Y)
                Y = Y/np.max(Y)
            VColors = np.array(np.round(255.0*cmap(Y)[:, 0:3]), dtype=np.int64)
            saveOffFileExternal("spectrum%i.off"%k, VPos, VColors, ITris)
    
    def doLowpass(self, evt):
        K = LOWPASS_K
        if (K >= self.mesh.VPos.shape[0]):
            K = self.mesh.VPos.shape[0]-1
        doLowpassFiltering(self.mesh, K)
        self.mesh.needsDisplayUpdate = True
        self.mesh.updateIndexDisplayList()
        self.Refresh()
    
    def doHeat(self, evt):
        K = HEAT_K
        if (K >= self.mesh.VPos.shape[0]):
            K = self.mesh.VPos.shape[0]-1
        (self.eigvalues, self.eigvectors) = getLaplacianSpectrum(self.mesh, K)
        self.GUIState = STATE_ANIMATEHEAT
        self.heatIdx = 0
        self.heat_ts = np.linspace(0, 1.0/self.eigvalues[1], 100)
        self.Refresh()

    def doHKS(self, evt):
        K = HKS_K
        if (K >= self.mesh.VPos.shape[0]):
            K = self.mesh.VPos.shape[0]-1
        heat = getHKS(self.mesh, K, HKS_T)
        heat = heat/np.max(heat)
        cmap = plt.get_cmap('jet')
        self.mesh.VColors = cmap(heat)[:, 0:3]
        self.mesh.needsDisplayUpdate = True
        self.Refresh()
    
    def doFlattening(self, evt):
        print self.laplacianSelections
        doFlattening(self.mesh, self.laplacianSelections)
        self.mesh.needsDisplayUpdate = True
        self.mesh.updateIndexDisplayList()
        self.Refresh()
    
    def doUVCoords(self, evt):
        (anchors, quadIdxs) = self.getAnchors()
        U = getTexCoords(self.mesh, quadIdxs)
        self.mesh.VTexCoords = U
        if self.mesh.texID == -1:
            self.mesh.texID = loadTexture("texture.png")
            print "Loaded texture: ", self.mesh.texID
        self.useTexture = True
        self.parent.useTextureCheckbox.SetValue(True)
        self.mesh.needsDisplayUpdate = True
        self.Refresh()
    
    ##Color callbacks
    def doSelectColorVertices(self, evt):
        if self.mesh:
            self.mesh.updateIndexDisplayList()
            self.GUIState = STATE_CHOOSECOLORVERTICES
            self.GUISubstate = COLORPICK_WAITING
            self.colorCurrentIdx = -1
            self.Refresh()
    
    def clearColorVertexSelection(self, evt):
        self.colorChoices.clear()
    
    def doInterpolateColors(self, evt):
        if not self.mesh:
            return
        (colors, colorsIdx) = self.getSelectedColors()
        if len(colorsIdx) == 0:
            return
        self.mesh.VColors = smoothColors(self.mesh, colors/255.0, colorsIdx)
        for k in range(3):
            self.mesh.VColors[self.mesh.VColors[:, k] < 0.0, k] = 0
            self.mesh.VColors[self.mesh.VColors[:, k] > 1.0 , k] = 1.0
        print self.mesh.VColors
        self.mesh.needsDisplayUpdate = True
        self.Refresh()
    
    def updateColorChoiceBuffers(self):
        #Remake color and vertex buffers
        if self.colorPosVBO:
            self.colorPosVBO.delete()
        if self.colorColorVBO:
            self.colorColorVBO.delete()
        N = len(self.colorChoices)
        Pos = np.zeros((N, 3))
        Colors = np.zeros((N, 3))
        i = 0
        for idx in self.colorChoices:
            Pos[i, :] = self.mesh.VPos[idx, :]
            Colors[i, :] = self.colorChoices[idx]/255.0
            i += 1
        self.colorPosVBO = vbo.VBO(np.array(Pos, dtype=np.float32))
        self.colorColorVBO = vbo.VBO(np.array(Colors, dtype=np.float32))
    
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
        elif self.GUIState == STATE_ANIMATEHEAT and self.eigvectors.size > 0 and self.eigvalues.size > 0:
            cmap = plt.get_cmap('jet')
            (anchors, anchorsIdx) = self.getAnchors()
            heat = getHeat(self.mesh, self.eigvalues, self.eigvectors, self.heat_ts[self.heatIdx], anchorsIdx)
#            if self.heatIdx == 0:
#                self.heatScale = np.max(heat)
#            heat = heat/self.heatScale
            self.mesh.VColors = cmap(heat)[:, 0:3]
            if self.mesh.VColorsVBO:
                self.mesh.VColorsVBO.delete()
            self.mesh.VColorsVBO = vbo.VBO(np.array(self.mesh.VColors, dtype=np.float32))
            self.drawMeshStandard()
            saveImageGL(self, "heat%i.png"%self.heatIdx)
            self.heatIdx += 1
            if self.heatIdx >= len(self.heat_ts):
                self.GUIState = STATE_NORMAL
            self.Refresh()
        elif self.GUIState == STATE_CHOOSELAPLACEVERTICES:
            if self.GUISubstate == CHOOSELAPLACE_WAITING:
                self.camera.gotoCameraFrame()
                glDisable(GL_LIGHTING)
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
                self.camera.gotoCameraFrame()
                self.mesh.renderGLIndices()
                pixel = glReadPixels(self.MousePos[0], self.MousePos[1], 1, 1, GL_RGBA, GL_UNSIGNED_BYTE)
                [R, G, B, A] = [int(pixel.encode("hex")[i*2:(i+1)*2], 16) for i in range(4)]
                idx = extractFromRGBA(R, G, B, 0) - 1
                if idx >= 0 and idx < len(self.mesh.vertices):
                    print idx
                    if idx in self.laplacianConstraints:
                        #De-select if it's already selected
                        self.laplaceCurrentIdx = -1
                        self.laplacianConstraints.pop(idx, None)
                        self.laplacianSelections.remove(idx)
                    else:
                        self.laplacianConstraints[idx] = np.array(self.mesh.VPos[self.mesh.vertices[idx].ID, :])
                        self.laplaceCurrentIdx = idx
                        self.laplacianSelections.append(idx)
                self.GUISubstate = CHOOSELAPLACE_WAITING
                self.Refresh()
                
        elif self.GUIState == STATE_CHOOSECOLORVERTICES:
            if self.GUISubstate == COLORPICK_WAITING:
                self.drawMeshStandard()
                if not self.colorPickTexID:
                    self.colorPickTexID = getColorPickingTexture()
                if self.colorPosVBO:
                    self.camera.gotoCameraFrame()
                    glDisable(GL_LIGHTING)
                    glEnableClientState(GL_VERTEX_ARRAY)
                    glEnableClientState(GL_COLOR_ARRAY)
                    self.colorPosVBO.bind()
                    glVertexPointerf(self.colorPosVBO)
                    self.colorColorVBO.bind()
                    glColorPointerf(self.colorColorVBO)
                    glPointSize(10)
                    glDrawArrays(GL_POINTS, 0, len(self.colorChoices))
                    self.colorPosVBO.unbind()
                    self.colorColorVBO.unbind()
                    glDisableClientState(GL_COLOR_ARRAY)
                    glDisableClientState(GL_VERTEX_ARRAY)
                drawColorPicker(self.size.width, self.size.height, self.colorPickTexID)
                
            elif self.GUISubstate == COLORPICK_PICKVERTEX:
                glClearColor(0.0, 0.0, 0.0, 0.0)
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
                self.camera.gotoCameraFrame()
                self.mesh.renderGLIndices()
                pixel = glReadPixels(self.MousePos[0], self.MousePos[1], 1, 1, GL_RGBA, GL_UNSIGNED_BYTE)
                [R, G, B, A] = [int(pixel.encode("hex")[i*2:(i+1)*2], 16) for i in range(4)]
                idx = extractFromRGBA(R, G, B, 0) - 1
                if idx >= 0 and idx < len(self.mesh.vertices):
                    if idx in self.colorChoices:
                        #De-select if it's already selected
                        self.colorCurrentIdx = -1
                        self.colorChoices.pop(idx, None)
                    else:
                        self.colorChoices[idx] = np.array([255, 255, 255])
                        self.colorCurrentIdx = idx
                    self.updateColorChoiceBuffers()
                self.GUISubstate = COLORPICK_WAITING
                self.Refresh()
            
            elif self.GUISubstate == COLORPICK_PICKCOLOR:
                self.drawMeshStandard()
                drawColorPicker(self.size.width, self.size.height, self.colorPickTexID)
                if self.colorCurrentIdx != -1 and self.MousePos[0] < 200 and self.MousePos[1] < 200:
                    pixel = glReadPixels(self.MousePos[0], self.MousePos[1], 1, 1, GL_RGBA, GL_UNSIGNED_BYTE)
                    [R, G, B, A] = [int(pixel.encode("hex")[i*2:(i+1)*2], 16) for i in range(4)]
                    self.colorChoices[self.colorCurrentIdx] = np.array([R, G, B])
                    self.updateColorChoiceBuffers()                    
                self.GUISubstate = COLORPICK_WAITING
                self.Refresh()
                
        self.SwapBuffers()

    def MouseDown(self, evt):
        state = wx.GetMouseState()
        if self.GUIState == STATE_CHOOSELAPLACEVERTICES:
            if state.ShiftDown():
                #Pick vertex for laplacian mesh constraints
                self.GUISubstate = CHOOSELAPLACE_PICKVERTEX
        elif self.GUIState == STATE_CHOOSECOLORVERTICES:
            if state.ShiftDown():
                self.GUISubstate = COLORPICK_PICKVERTEX
            elif state.CmdDown() or self.pressingC:
                self.GUISubstate = COLORPICK_PICKCOLOR
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
            if self.GUIState == STATE_CHOOSELAPLACEVERTICES and (state.CmdDown() or self.pressingC) and self.laplaceCurrentIdx in self.laplacianConstraints:
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
    (ID_LOADDATASET, ID_SAVEDATASET, ID_SAVEDATASETMETERS, ID_SAVESCREENSHOT, ID_CONNECTEDCOMPONENTS, ID_SPLITFACES, ID_TRUNCATE, ID_FILLHOLES, ID_GEODESICDISTANCES, ID_PRST, ID_INTERPOLATECOLORS, ID_SAVEROTATINGSCREENSOTS, ID_SAVELIGHTINGSCREENSHOTS, ID_SELECTLAPLACEVERTICES, ID_CLEARLAPLACEVERTICES, ID_SOLVEWITHCONSTRAINTS, ID_MEMBRANEWITHCONSTRAINTS, ID_GETHKS, ID_GETHEATFLOW, ID_LAPLACIANSMOOTH, ID_LAPLACIANSHARPEN, ID_MINIMALSURFACE, ID_ESTIMATECURVATURE, ID_GETSPECTRUM, ID_DOLOWPASS, ID_DOHEAT, ID_DOHKS, ID_SELECTCOLORVERTICES, ID_CLEARCOLORVERTICES, ID_INTERPOLATECOLORS, ID_DOFLATTENING, ID_DOUVCOORDS) = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32)
    
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
        
#        menuEstimateCurvature = laplacianMenu.Append(MeshViewerFrame.ID_ESTIMATECURVATURE, "&Estimate Mean Curvature", "Estimate Mean Curvature")
#        self.Bind(wx.EVT_MENU, self.glcanvas.estimateCurvature, menuEstimateCurvature)
        
        menuLaplacianSmooth = laplacianMenu.Append(MeshViewerFrame.ID_LAPLACIANSMOOTH, "&Laplacian Smooth", "Laplacian Smooth")
        self.Bind(wx.EVT_MENU, self.glcanvas.doLaplacianSmooth, menuLaplacianSmooth)
        
        menuLaplacianSharpen = laplacianMenu.Append(MeshViewerFrame.ID_LAPLACIANSHARPEN, "&Laplacian Sharpen", "Laplacian Sharpen")
        self.Bind(wx.EVT_MENU, self.glcanvas.doLaplacianSharpen, menuLaplacianSharpen)
        
        menuMinimalSurface = laplacianMenu.Append(MeshViewerFrame.ID_MINIMALSURFACE, "&Minimal Surface", "Minimal Surface")
        self.Bind(wx.EVT_MENU, self.glcanvas.doMinimalSurface, menuMinimalSurface)
        
        menuGetSpectrum = laplacianMenu.Append(MeshViewerFrame.ID_GETSPECTRUM, "&Get Spectrum", "Get Spectrum")
        self.Bind(wx.EVT_MENU, self.glcanvas.getSpectrum, menuGetSpectrum)

        menuDoLowpass = laplacianMenu.Append(MeshViewerFrame.ID_DOLOWPASS, "&Do Lowpass", "Do Lowpass")
        self.Bind(wx.EVT_MENU, self.glcanvas.doLowpass, menuDoLowpass)
        
        menuDoHeat = laplacianMenu.Append(MeshViewerFrame.ID_DOHEAT, "&Do Heat Flow Simulation", "Do Heat Flow Simulation")
        self.Bind(wx.EVT_MENU, self.glcanvas.doHeat, menuDoHeat)

        menuDoHKS = laplacianMenu.Append(MeshViewerFrame.ID_DOHKS, "&Compute Heat Kernel Signature", "Do Heat Heat Kernel Signature")
        self.Bind(wx.EVT_MENU, self.glcanvas.doHKS, menuDoHKS)
        
        menuDoFlattening = laplacianMenu.Append(MeshViewerFrame.ID_DOFLATTENING, "&Do Flattening", "Do Flattening")
        self.Bind(wx.EVT_MENU, self.glcanvas.doFlattening, menuDoFlattening)
        
        menuDoUVCoords = laplacianMenu.Append(MeshViewerFrame.ID_DOUVCOORDS, "&Compute UV Coordinates", "Compute UV Coordinates")
        self.Bind(wx.EVT_MENU, self.glcanvas.doUVCoords, menuDoUVCoords) 
        
        #####Color Selection Menu
        colorMenu = wx.Menu()
        menuSelectColorVertices = colorMenu.Append(MeshViewerFrame.ID_SELECTCOLORVERTICES, "&Select Color Vertices", "Select Color Vertices")
        self.Bind(wx.EVT_MENU, self.glcanvas.doSelectColorVertices, menuSelectColorVertices)
        
        menuClearColorVertices = colorMenu.Append(MeshViewerFrame.ID_CLEARCOLORVERTICES, "&Clear color vertex selection", "Clear Color Vertex Selection")
        self.Bind(wx.EVT_MENU, self.glcanvas.clearColorVertexSelection, menuClearColorVertices)
        
        menuInterpolateColors = colorMenu.Append(MeshViewerFrame.ID_INTERPOLATECOLORS, "&Interpolate Colors", "Interpolate Colors")
        self.Bind(wx.EVT_MENU, self.glcanvas.doInterpolateColors, menuInterpolateColors)
        
        # Creating the menubar.
        menuBar = wx.MenuBar()
        menuBar.Append(filemenu,"&File") # Adding the "filemenu" to the MenuBar
        menuBar.Append(laplacianMenu,"&MeshLaplacian") # Adding the "filemenu" to the MenuBar
        menuBar.Append(colorMenu,"&MeshColoring") #Adding the coloring menu
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
        dlg = wx.FileDialog(self, "Choose a file", ".", "", "OFF files (*.off)|*.off|TOFF files (*.toff)|*.toff|OBJ files (*.obj)|*.obj", wx.FD_OPEN)
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
            self.glcanvas.clearAllSelections()
            self.glcanvas.Refresh()
        dlg.Destroy()
        return

    def OnSaveMesh(self, evt):
        dlg = wx.FileDialog(self, "Choose a file", ".", "", "*", wx.FD_SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetFilename()
            dirname = dlg.GetDirectory()
            filepath = os.path.join(dirname, filename)
            self.glcanvas.mesh.saveFile(filepath, True)
            self.glcanvas.Refresh()
        dlg.Destroy()
        return     
        
    def OnSaveScreenshot(self, evt):
        dlg = wx.FileDialog(self, "Choose a file", ".", "", "*", wx.FD_SAVE)
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