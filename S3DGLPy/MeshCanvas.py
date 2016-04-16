#Contains a basic class for viewing a mesh, wrapping around glcanvas.  Includes
#Basic functionality for setting up lighting/cameras and handling mouse events
#for changing viewing parameters
from OpenGL.GL import *
from OpenGL.arrays import vbo
import wx
from wx import glcanvas

from Primitives3D import *
from PolyMesh import *
from LaplacianMesh import *
from Cameras3D import *

DEFAULT_SIZE = wx.Size(1200, 800)
DEFAULT_POS = wx.Point(10, 10)

def saveImageGL(mvcanvas, filename):
    view = glGetIntegerv(GL_VIEWPORT)
    img = wx.EmptyImage(view[2], view[3] )
    pixels = glReadPixels(0, 0, view[2], view[3], GL_RGB,
                     GL_UNSIGNED_BYTE)
    img.SetData( pixels )
    img = img.Mirror(False)
    img.SaveFile(filename, wx.BITMAP_TYPE_PNG)

def saveImage(canvas, filename):
    s = wx.ScreenDC()
    w, h = canvas.size.Get()
    b = wx.EmptyBitmap(w, h)
    m = wx.MemoryDCFromDC(s)
    m.SelectObject(b)
    m.Blit(0, 0, w, h, s, 70, 0)
    m.SelectObject(wx.NullBitmap)
    b.SaveFile(filename, wx.BITMAP_TYPE_PNG)

class BasicMeshCanvas(glcanvas.GLCanvas):
    def __init__(self, parent):
        attribs = (glcanvas.WX_GL_RGBA, glcanvas.WX_GL_DOUBLEBUFFER, glcanvas.WX_GL_DEPTH_SIZE, 24)
        glcanvas.GLCanvas.__init__(self, parent, -1, attribList = attribs)    
        self.context = glcanvas.GLContext(self)
        
        self.parent = parent
        #Camera state variables
        self.size = self.GetClientSize()
        #self.camera = MouseSphericalCamera(self.size.x, self.size.y)
        self.camera = MousePolarCamera(self.size.width, self.size.height)
        
        #Main state variables
        self.MousePos = [0, 0]
        self.bbox = BBox3D()  
        
        #Face mesh variables and manipulation variables
        self.mesh = None
        self.meshCentroid = None
        self.displayMeshFaces = True
        self.displayMeshEdges = False
        self.displayMeshVertices = True
        self.displayVertexNormals = False
        self.displayFaceNormals = False
        self.useLighting = True
        self.useTexture = False
        
        self.GLinitialized = False
        #GL-related events
        wx.EVT_ERASE_BACKGROUND(self, self.processEraseBackgroundEvent)
        wx.EVT_SIZE(self, self.processSizeEvent)
        wx.EVT_PAINT(self, self.processPaintEvent)
        #Mouse Events
        wx.EVT_LEFT_DOWN(self, self.MouseDown)
        wx.EVT_LEFT_UP(self, self.MouseUp)
        wx.EVT_RIGHT_DOWN(self, self.MouseDown)
        wx.EVT_RIGHT_UP(self, self.MouseUp)
        wx.EVT_MIDDLE_DOWN(self, self.MouseDown)
        wx.EVT_MIDDLE_UP(self, self.MouseUp)
        wx.EVT_MOTION(self, self.MouseMotion)
    
    def initMeshBBox(self):
        if self.mesh:
            self.bbox = self.mesh.getBBox()
            print "Mesh BBox: %s\n"%self.bbox
            self.camera.centerOnBBox(self.bbox, theta = -math.pi/2, phi = math.pi/2)
        
    def viewFromFront(self, evt):
        self.camera.centerOnBBox(self.bbox, theta = -math.pi/2, phi = math.pi/2)
        self.Refresh()
    
    def viewFromTop(self, evt):
        self.camera.centerOnBBox(self.bbox, theta = -math.pi/2, phi = 0)
        self.Refresh()
    
    def viewFromSide(self, evt):
        self.camera.centerOnBBox(self.bbox, theta = -math.pi, phi = math.pi/2)
        self.Refresh()
    
    def processEraseBackgroundEvent(self, event): pass #avoid flashing on MSW.
    
    def processSizeEvent(self, event):
        self.size = self.GetClientSize()
        glViewport(0, 0, self.size.width, self.size.height)
        #Update camera parameters based on new size
        self.camera = MousePolarCamera(self.size.width, self.size.height)
        self.camera.centerOnBBox(self.bbox, math.pi/2, math.pi/2)

    def processPaintEvent(self, event):
        dc = wx.PaintDC(self)
        self.SetCurrent(self.context)
        if not self.GLinitialized:
            self.initGL()
            self.GLinitialized = True
        self.repaint()

    def drawMeshStandard(self):
        glEnable(GL_LIGHTING)
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, [0.8, 0.8, 0.8, 1.0])
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, [0.2, 0.2, 0.2, 1.0])
        glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, 64)
        #Set up modelview matrix
        self.camera.gotoCameraFrame()
        glLightfv(GL_LIGHT0, GL_POSITION, np.array([0, 0, 0, 1]))
        self.mesh.renderGL(self.displayMeshEdges, self.displayMeshVertices, self.displayMeshFaces, self.displayVertexNormals, self.displayFaceNormals, self.useLighting, self.useTexture)
    
    def setupPerspectiveMatrix(self, nearDist = -1, farDist = -1):
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        if nearDist == -1:
            farDist = self.camera.eye - self.bbox.getCenter()
            farDist = np.sqrt(farDist.dot(farDist)) + self.bbox.getDiagLength()
            nearDist = farDist/500.0
        gluPerspective(180.0*self.camera.yfov/M_PI, float(self.size.x)/self.size.y, nearDist, farDist)
    
    def repaint(self):
        self.setupPerspectiveMatrix()
        glClearColor(0.0, 0.0, 0.0, 0.0)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        
        if self.mesh:
            self.drawMeshStandard()
        self.SwapBuffers()
    
    def initGL(self):        
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, [0.2, 0.2, 0.2, 1.0])
        glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE)
        glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
        glEnable(GL_LIGHT0)
        glLightfv(GL_LIGHT1, GL_DIFFUSE, [0.5, 0.5, 0.5, 1.0])
        glEnable(GL_LIGHT1)
        glEnable(GL_NORMALIZE)
        glEnable(GL_LIGHTING)
        glEnable(GL_DEPTH_TEST)

    def handleMouseStuff(self, x, y):
        #Invert y from what the window manager says
        y = self.size.height - y
        self.MousePos = [x, y]

    def MouseDown(self, evt):
        state = wx.GetMouseState()
        x, y = evt.GetPosition()
        self.CaptureMouse()
        self.handleMouseStuff(x, y)
        self.Refresh()
    
    def MouseUp(self, evt):
        x, y = evt.GetPosition()
        self.handleMouseStuff(x, y)
        self.ReleaseMouse()
        self.Refresh()

    def MouseMotion(self, evt):
        state = wx.GetMouseState()
        x, y = evt.GetPosition()
        [lastX, lastY] = self.MousePos
        self.handleMouseStuff(x, y)
        dX = self.MousePos[0] - lastX
        dY = self.MousePos[1] - lastY
        if evt.Dragging():
            #Translate/rotate shape
            if evt.MiddleIsDown():
                self.camera.translate(dX, dY)
            elif evt.RightIsDown():
                self.camera.zoom(-dY)#Want to zoom in as the mouse goes up
            elif evt.LeftIsDown():
                self.camera.orbitLeftRight(dX)
                self.camera.orbitUpDown(dY)
        self.Refresh() 

if __name__ == '__main__':
    app = wx.PySimpleApp()
    frame = wx.Frame(None, wx.ID_ANY, "Basic Mesh Canvas", DEFAULT_POS, DEFAULT_SIZE)
    g = BasicMeshCanvas(frame)
    g.mesh = getDodecahedronMesh()
    g.initMeshBBox()
    frame.canvas = g
    frame.Show()
    app.MainLoop()
    app.Destroy()
