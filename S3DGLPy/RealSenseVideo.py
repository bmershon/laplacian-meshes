import matplotlib.pyplot as plt
import matplotlib.image as mpimg
#from mpl_toolkits.mplot3d import Axes3D
import scipy.io as sio
import scipy.misc
from scipy import ndimage
import numpy as np
import numpy.linalg as linalg
from PolyMesh import *
from LaplacianMesh import *
from OpenGL.arrays import vbo
from OpenGL.GL import *
from VideoTools import *
import time

def imreadf(filename):
    #Read in file, converting image byte array to little endian float
    I = scipy.misc.imread(filename)
    #Image is stored in BGRA format so convert to RGBA
    I = I[:, :, [2, 1, 0, 3]]
    shape = I.shape
    I = I.flatten()
    IA = bytearray(I.tolist())
    I = np.fromstring(IA.__str__(), dtype=np.dtype('<f4'))
    return np.reshape(I, shape[0:2])  

def imwritef(I, filename):
    IA = I.flatten().tolist()
    IA = struct.pack("%if"%len(IA), *IA)
    IA = np.fromstring(IA, dtype=np.uint8)
    IA = IA.reshape([I.shape[0], I.shape[1], 4]) ##Tricky!!  Numpy is "low-order major" and the order I have things in is 4bytes per pixel, then columns, then rows.  These are specified in reverse order
    print "IA.shape = ", IA.shape
    #Convert from RGBA format to BGRA format like the real sense saver did
    IA = IA[:, :, [2, 1, 0, 3]]
    scipy.misc.imsave(filename, IA)

#Use the uv coorinates to map into the array of colors "C"
#using bilinear interpolation.  Out of bounds values are by default gray
def getColorsFromMap(u, v, C, mask):
    #The uv coordinates are in [0, 1]x[0, 1] so scale to [0, height]x[0, width]
    [H, W] = [C.shape[0], C.shape[1]]
    thisu = u[mask > 0]*W
    thisv = v[mask > 0]*H
    #Round out of bounds indices to the edge for now
    thisu[thisu >= W] = W-1
    thisu[thisu < 0] = 0
    thisv[thisv >= H] = H-1
    thisv[thisv < 0] = 0
    N = len(thisu)
    #Do bilinear interpolation on grid
    m = mask.flatten()
    CAll = np.reshape(C, [C.shape[0]*C.shape[1], C.shape[2]])
    c = CAll[m > 0, :]
    
    #utop, ubottom
    ul = np.array(np.floor(thisu), dtype=np.int64)
    ur = ul + 1
    ur[ur >= W] = W-1
    #vteft, vbight
    vt = np.array(np.floor(thisv), dtype=np.int64)
    vb = vt + 1
    vb[vb >= H] = H-1
    c = (ur-thisu)[:, None]*(CAll[vt*W+ul, :]*(vb-thisv)[:, None] + CAll[vb*W+ul, :]*(thisv-vt)[:, None]) + (thisu-ul)[:, None]*(CAll[vt*W+ur, :]*(vb-thisv)[:, None] + CAll[vb*W+ur, :]*(thisv-vt)[:, None])
    #Set out of bounds pixels to gray (since color/depth aren't aligned perfectly)
    thisu = u[mask > 0]*W
    thisv = v[mask > 0]*H
    c[thisu >= W, :] = np.array([0.5, 0.5, 0.5])
    c[thisu < 0, :] = np.array([0.5, 0.5, 0.5])
    c[thisv >= H, :] = np.array([0.5, 0.5, 0.5])
    c[thisv < 0, :] = np.array([0.5, 0.5, 0.5])
    return c
    

def getFrame(foldername, index, loadColor = True, plotFrame = False):
    depthFile = "%s/B-depth-float%i.png"%(foldername, index)
    xyFile = "%s/B-cloud%i.png"%(foldername, index)
    Z = imreadf(depthFile)
    XYZ = imreadf(xyFile)
    X = XYZ[:, 0:-1:3]
    Y = XYZ[:, 1:-1:3]
    uvname = "%s/B-depth-uv%i.png"%(foldername, index)
    u = np.zeros(Z.shape)
    v = np.zeros(Z.shape)
    C = 0.5*np.ones((Z.shape[0], Z.shape[1], 3)) #Default gray
    loadedColor = False
    if loadColor and os.path.exists(uvname):
        uv = imreadf(uvname)
        u = uv[:, 0::2]
        v = uv[:, 1::2]
        C = scipy.misc.imread("%s/B-color%i.png"%(foldername, index)) / 255.0
        loadedColor = True
    if plotFrame:
        x = X[Z > 0]
        y = Y[Z > 0]
        z = Z[Z > 0]
        c = getColorsFromMap(u, v, C, Z > 0)
        fig = plt.figure()
        #ax = Axes3D(fig)
        plt.scatter(x, y, 30, c)
        plt.show()
    return [X, Y, Z, C, u, v, loadedColor]

class RealSenseVideo(object):
    def __init__(self):
        self.mesh = None
        self.VPosVBOs = []
        self.VNormalsVBOs = []
        self.Xs = np.zeros((0, 0, 0))
        self.Ys = np.zeros((0, 0, 0))
        self.Zs = np.zeros((0, 0, 0))
        self.Cs = np.zeros((0, 0, 0, 0)) #Color frames
        self.us = np.zeros((0, 0, 0)) #depth to color map horiz coordinate
        self.vs = np.zeros((0, 0, 0)) #depth to color map vert coordinate
        self.rawAmplification = False
        self.loadedColor = False
        
        #Amplification variables
        self.Ls = []
        self.M_invs = []
        self.solvers = []
        self.gamma = 1.0
        self.origDeltaCoords = np.zeros((0, 0))
        self.ampDeltaCoords = np.zeros((0, 0))

    def loadVideo(self, foldername, NFrames, loadColor = True):
        if NFrames <= 0:
            return
        shape = getFrame(foldername, 0)[0].shape
        Xs = np.zeros((shape[0], shape[1], NFrames))
        Ys = np.zeros(Xs.shape)
        Zs = np.zeros(Xs.shape)
        Cs = 0.5*np.ones([Xs.shape[0], Xs.shape[1], 3, NFrames])
        us = np.zeros(Xs.shape)
        vs = np.zeros(Xs.shape)
        for i in range(NFrames):
            print "Loading %s frame %i"%(foldername, i)
            [Xs[:, :, i], Ys[:, :, i], Zs[:, :, i], C, us[:, :, i], vs[:, :, i], self.loadedColor] = getFrame(foldername, i, loadColor)
            if self.loadedColor:
                if i == 0:
                    Cs = np.zeros((C.shape[0], C.shape[1], C.shape[2], NFrames))
                Cs[:, :, :, i] = C
        self.Xs = Xs
        self.Ys = Ys
        self.Zs = Zs
        self.Cs = Cs
        self.us = us
        self.vs = vs

    def makeMask(self):
        #Narrow down to pixels which measured something every frame
        counts = np.sum(self.Zs > 0, 2)
        self.Mask = (counts == self.Zs.shape[2])
        #Find biggest connected component out of remaining pixels
        ILabel, NLabels = ndimage.label(self.Mask)
        idx = np.argmax(ndimage.sum(self.Mask, ILabel, range(NLabels+1)))
        self.Mask = (ILabel == idx)
        plt.imshow(self.Mask)
        plt.show()

    #Actually sets up all of the vertex, face, and adjacency structures in the 
    #PolyMeshObject
    def makeMeshSlow(self):
        if self.Xs.shape[2] == 0:
            return
        X = self.Xs[:, :, 0]
        Y = self.Ys[:, :, 0]
        Z = self.Zs[:, :, 0]
        mesh = PolyMesh()
        Vs = [[None]*self.Mask.shape[1] for i in range(self.Mask.shape[0])]
        #Add vertices
        for i in range(self.Mask.shape[0]):
            for j in range(self.Mask.shape[1]):
                if self.Mask[i, j]:
                    Vs[i][j] = mesh.addVertex(np.array([X[i, j], Y[i, j], -Z[i, j]]))
        #Add triangles on grid
        for i in range(self.Mask.shape[0]-1):
            print i
            for j in range(self.Mask.shape[1]-1):
                if self.Mask[i, j] and self.Mask[i+1, j] and self.Mask[i, j+1]:
                    mesh.addFace([Vs[i][j], Vs[i+1][j], Vs[i][j+1]])
                if self.Mask[i+1][j] and self.Mask[i+1][j+1] and self.Mask[i][j+1]:
                    mesh.addFace([Vs[i+1][j], Vs[i+1][j+1], Vs[i][j+1]])
        mesh.updateTris()
        self.mesh = mesh
    
    #Bypasses a lot of the PolyMesh structure to fill in vertex, triangle, and
    #normal buffer only, but the mesh object has to be treated with care when
    #used since the structure is not there
    def makeMesh(self):
        if self.Xs.shape[2] == 0:
            return
        X = self.Xs[:, :, 0]
        Y = self.Ys[:, :, 0]
        Z = self.Zs[:, :, 0]
        mesh = PolyMesh()
        #Come up with vertex indices in the mask
        Mask = np.array(self.Mask, dtype=np.int32)
        nV = np.sum(Mask)
        Mask[self.Mask > 0] = np.arange(nV) + 1
        Mask = Mask - 1
        VPos = np.zeros((nV, 3))
        VPos[:, 0] = X[self.Mask > 0]
        VPos[:, 1] = Y[self.Mask > 0]
        VPos[:, 2] = -Z[self.Mask > 0]
        #Add lower right triangle
        v1 = Mask[0:-1, 0:-1].flatten()
        v2 = Mask[1:, 0:-1].flatten()
        v3 = Mask[1:, 1:].flatten()
        N = v1.size
        ITris1 = np.concatenate((np.reshape(v1, [N, 1]), np.reshape(v2, [N, 1]), np.reshape(v3, [N, 1])), 1)
        #Add upper right triangle
        v1 = Mask[0:-1, 0:-1].flatten()
        v2 = Mask[1:, 1:].flatten()
        v3 = Mask[0:-1, 1:].flatten()
        N = v1.size
        ITris2 = np.concatenate((np.reshape(v1, [N, 1]), np.reshape(v2, [N, 1]), np.reshape(v3, [N, 1])), 1)
        ITris = np.concatenate((ITris1, ITris2), 0)
        #Only retain triangles which have all three points
        ITris = ITris[np.sum(ITris == -1, 1) == 0, :]
        mesh.VPos = VPos
        mesh.ITris = ITris
        mesh.VColors = 0.5*np.ones(mesh.VPos.shape)
        mesh.updateNormalBuffer()
        mesh.VPosVBO = vbo.VBO(np.array(mesh.VPos, dtype=np.float32))
        mesh.VNormalsVBO = vbo.VBO(np.array(mesh.VNormals, dtype=np.float32))
        mesh.VColorsVBO = vbo.VBO(np.array(mesh.VColors, dtype=np.float32))
        mesh.IndexVBO = vbo.VBO(mesh.ITris, target=GL_ELEMENT_ARRAY_BUFFER)
        mesh.needsDisplayUpdate = False
        self.mesh = mesh

    #Compute the delta coordinates and solvers for all frames of the video
    def computeLaplacians(self, gamma):
        if not self.mesh:
            self.makeMesh()
        self.gamma = gamma
        self.solvers = []
        self.Ls = []
        self.M_invs = []
        NFrames = self.Xs.shape[2]
        self.origDeltaCoords = np.zeros((self.mesh.VPos.shape[0]*6, NFrames))
        for i in range(NFrames):
            print "Getting laplace embedding frame %i of %i\n"%(i, NFrames)
            VPos = np.zeros(self.mesh.VPos.shape)
            [X, Y, Z] = [self.Xs[:, :, i], self.Ys[:, :, i], self.Zs[:, :, i]]
            VPos[:, 0] = X[self.Mask > 0]
            VPos[:, 1] = Y[self.Mask > 0]
            VPos[:, 2] = -Z[self.Mask > 0]
            anchorsIdx = np.arange(VPos.shape[0])
            (L, M_inv, solver, deltaCoords) = makeLaplacianMatrixSolverIGLSoft(VPos, self.mesh.ITris, anchorsIdx, gamma, makeSolver = False)
            deltaCoords = np.concatenate((deltaCoords, VPos), 0)
            self.origDeltaCoords[:, i] = deltaCoords.flatten()
    
    #Do amplification, launching GUI to help choose singular vectors
    def doAmplification(self, W, alpha):
        self.rawAmplification = False
#        #Step 1: Get the right singular vectors
#        Mu = tde_mean(self.origDeltaCoords, W)
#        (Y, S) = tde_rightsvd(self.origDeltaCoords, W, Mu)
#        
#        #Step 2: Choose which components to amplify by a visual inspection
#        chooser = PCChooser(Y, (alpha, DEFAULT_NPCs))
#        Alpha = chooser.getAlphaVec(alpha)
        
        #Step 3: Perform the amplification
        #Add the amplified delta coordinates back to the original delta coordinates
        self.ampDeltaCoords = self.origDeltaCoords + subspace_tde_amplification(self.origDeltaCoords, self.origDeltaCoords.shape, W, alpha*np.ones((1,DEFAULT_NPCs)), self.origDeltaCoords.shape[1]/2)
        
#        self.ampDeltaCoords = self.origDeltaCoords + tde_amplifyPCs(self.origDeltaCoords, W, Mu, Y, Alpha)
#        print 'normalized error:',(linalg.norm(self.ampDeltaCoords-other_delta_coords)/linalg.norm(self.ampDeltaCoords))
#        print "Finished Amplifying"

    #Perform an amplification on the raw XYZ coordinates
    def doAmplificationRaw(self, W, alpha):
        self.rawAmplification = True
        NFrames = self.Xs.shape[2]
        origPos = np.zeros((self.mesh.VPos.shape[0]*3, NFrames))
        for i in range(NFrames):
            VPos = np.zeros(self.mesh.VPos.shape)
            [X, Y, Z] = [self.Xs[:, :, i], self.Ys[:, :, i], self.Zs[:, :, i]]
            VPos[:, 0] = X[self.Mask > 0]
            VPos[:, 1] = Y[self.Mask > 0]
            VPos[:, 2] = -Z[self.Mask > 0]
            origPos[:, i] = VPos.flatten()
        self.newPos = origPos + subspace_tde_amplification(origPos, origPos.shape, W, alpha, origPos.shape[1])

    #Copy in vertices and compute vertex normals for a frame
    #Used to update visualization for original video and amplified video
    def swapInFrame(self, i, playOrigVideo = True, updateNormals = True):
        if not self.mesh:
            self.makeMesh()
        if i > self.Xs.shape[2]:
            print "Warning: Trying to load in frame %i that is beyond range of %i frames"%(i, self.Xs.shape[2])
            return
        m = self.mesh
        VPos = np.zeros(m.VPos.shape)
        Mask = self.Mask
        [X, Y, Z] = [self.Xs[:, :, i], self.Ys[:, :, i], self.Zs[:, :, i]]
        VPos[:, 0] = X[Mask > 0]
        VPos[:, 1] = Y[Mask > 0]
        VPos[:, 2] = -Z[Mask > 0]
        if playOrigVideo:
            m.VPos = VPos
        elif self.rawAmplification:
            m.VPos = np.reshape(self.newPos[:, i], VPos.shape)
        else:
            deltaCoordsAmp = np.reshape(self.ampDeltaCoords[:, i], [VPos.shape[0]*2, 3])
            #Run the solver on the amplified delta coordinates and determine the frame vertices
            (L, M_inv, solver, deltaCoords) = makeLaplacianMatrixSolverIGLSoft(VPos, m.ITris, np.arange(VPos.shape[0]), self.gamma, makeSolver = True)
            deltaCoords = deltaCoordsAmp[0:VPos.shape[0], :]
            posAmp = deltaCoordsAmp[VPos.shape[0]:, :]
            m.VPos = solveLaplacianMatrixIGLSoft(solver, L, M_inv, deltaCoordsAmp, np.arange(VPos.shape[0]), VPos, self.gamma)      
        
        if m.VPosVBO:
            m.VPosVBO.delete()
        if updateNormals:
            m.updateNormalBuffer()
            if m.VNormalsVBO:
                m.VNormalsVBO.delete()
            m.VPosVBO = vbo.VBO(np.array(m.VPos, dtype=np.float32))
            m.VNormalsVBO = vbo.VBO(np.array(m.VNormals, dtype=np.float32))
        if self.loadedColor:
            if m.VColorsVBO:
                m.VColorsVBO.delete()
            c = getColorsFromMap(self.us[:, :, i], self.vs[:, :, i], self.Cs[:, :, :, i], self.Mask)
            m.VColors = c
            print "Getting colors"
            m.VColorsVBO = vbo.VBO(np.array(c, dtype=np.float32))
    
    def saveFrame(self, i, filename, playOrigVideo = True):
        self.swapInFrame(i, playOrigVideo)
        m = self.mesh
        saveOffFileExternal(filename, m.VPos, m.VColors, m.ITris)

#Comparing real sense to kinect frames
if __name__ == '__main__':
    v1 = RealSenseVideo()
    v1.loadVideo("3DVideos/Chris_Neck", 1)
    v1.Mask = v1.Zs[:, :, 0] > 0
    v1.saveFrame(0, "RealSense.off")
    
    v2 = RealSenseVideo()
    v2.loadVideo("3DVideos/Chris_Head_Color_KinectV2", 1)
    v2.Mask = v2.Zs[:, :, 0] > 0
    v2.saveFrame(0, "KinectV2.off")

#Testing colored point cloud with UV map
if __name__ == '__main__2':
    [X, Y, Z, C, u, v, loadedColor] = getFrame("Chris_Neck_Color", 0, True)
    x = X[Z > 0]
    y = Y[Z > 0]
    z = Z[Z > 0]
    tic = time.time()
    c = getColorsFromMap(u, v, C, Z > 0)
    toc = time.time()
    print "Grid interp elapsed time ", toc - tic
    VPos = np.array([x.flatten(), y.flatten(), z.flatten()]).T
    saveOffFileExternal("ColorPC.off", VPos, c, np.zeros((0, 3)))

