#Programmer: Chris Tralie
#Purpose: To provide a simple polygon mesh class on top of numpy and PyOpenGL
#capable of reading and writing off files, as well as rendering meshes and
#performing some simple geometric and topological manipulations
from Primitives3D import *
from OpenGL.GL import *
from OpenGL.arrays import vbo
import OpenGL.GL.shaders
import sys
import re
import numpy as np
import numpy.linalg as linalg
import struct

POINT_SIZE = 7

def loadTexture(filename):
    try:
        from PIL.Image import open as imgopen
    except ImportError, err:
        from Image import open as imgopen
    im = imgopen(filename)
    try:
        im = im.convert('RGB')
        ix, iy, image = im.size[0], im.size[1], im.tobytes("raw", "RGBA", 0, -1)
    except SystemError:
        ix, iy, image = im.size[0], im.size[1], im.tobytes("raw", "RGBX", 0, -1)
    assert ix*iy*4 == len(image), """Unpacked image size for texture is incorrect"""
    
    texID = glGenTextures(1)
    glBindTexture(GL_TEXTURE_2D, texID)
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
    glTexImage2D(GL_TEXTURE_2D, 0, 3, ix, iy, 0, GL_RGBA, GL_UNSIGNED_BYTE, image)
    return texID

#Used to help sort edges in CCW order
class EdgesCCWComparator(object):
    def __init__(self, VCenter, N):
        self.VCenter = VCenter #Common point of all edges
        self.N = N
        self.mesh = VCenter.mesh
    
    def compare(self, e1, e2):
        V1 = e1.vertexAcross(self.VCenter)
        V2 = e2.vertexAcross(self.VCenter)
        a = V1.getPos() - VCenter.getPos()
        b = V2.getPos() - VCenter.getPos()
        triNormal = np.cross(a, b)
        dot = triNormal.dot(self.N)
        if dot > 0:
            return 1
        elif dot == 0:
            return 0
        return -1

class MeshVertex(object):
    def __init__(self, mesh, pos, color, texCoords, ID):
        self.mesh = mesh
        self.ID = ID
        #Set color and position in buffers
        if self.mesh.VPos.shape[0] <= ID:
            #If the client side vertex array needs to be expanded
            NPad = ID - self.mesh.VPos.shape[0] + 1
            self.mesh.VPos = np.concatenate((self.mesh.VPos, np.zeros((NPad, 3))), 0)
            self.mesh.VColors = np.concatenate((self.mesh.VColors, np.zeros((NPad, 3))), 0)
            self.mesh.VTexCoords = np.concatenate((self.mesh.VTexCoords, np.zeros((NPad, 2))), 0)
        self.mesh.VPos[ID, :] = pos
        self.mesh.VColors[ID, :] = color
        self.mesh.VTexCoords[ID, :] = texCoords
        self.edges = set() #Store reference to all emanating edges
        #NOTE: Edges are not guaranteed to be in any particular order
        self.component = -1 #Which connected component it's in
    
    def getPos(self):
        return self.mesh.VPos[self.ID, :]
    
    def getVertexNeighbors(self):
        ret = [0]*len(self.edges)
        i = 0
        for edge in self.edges:
            ret[i] = edge.vertexAcross(self)
            i = i+1
        return ret
    
    #Return a set of all faces attached to this vertex
    def getAttachedFaces(self):
        ret = set()
        i = 0
        for edge in self.edges:
            if (edge.f1 != None):
                ret.add(edge.f1)
            if (edge.f2 != None):
                ret.add(edge.f2)
        return ret
    
    #Return the area of the one-ring faces attached to this vertex
    def getOneRingArea(self):
        faces = self.getAttachedFaces()
        ret = 0.0
        for f in faces:
            ret = ret + f.getArea()
        return ret
    
    #Get an estimate of the vertex normal by taking a weighted
    #average of normals of attached faces
    def getNormal(self):
        faces = self.getAttachedFaces()
        totalArea = 0.0;
        normal = np.array([0, 0, 0])
        for f in faces:
            w = f.getArea()
            totalArea = totalArea + w
            normal = normal + w*f.getNormal()
        if totalArea == 0:
            return normal
        return (1.0/totalArea)*normal
    
    #Sort the edges so that they are in CCW order when projected onto
    #the plane defined by the vertex's normal
    def getCCWSortedEdges(self):
        comparator = EdgesCCWComparator(self, self.getNormal())
        return sorted(self.edges, cmp = comparator.compare)

class MeshFace(object):
    def __init__(self, mesh, ID):
        self.mesh = mesh
        self.ID = ID
        self.edges = [] #Store edges in CCW order
        self.startV = 0 #Vertex that starts it off

    def flipNormal(self):
        #Reverse the specification of the edges to make the normal
        #point in the opposite direction
        self.edges.reverse()
        self.normal = np.zeros((0, 0)) #Invalidate current cached normal

    #Return a list of vertices on the face in CCW order
    def getVertices(self):
        ret = [0]*len(self.edges)
        v = self.startV
        for i in range(0, len(self.edges)):
            ret[i] = v
            v = self.edges[i].vertexAcross(v)
        return ret
    
    #Return the vertices' positions in an N x 3 matrix
    def getVerticesPos(self):
        return self.mesh.VPos[[v.ID for v in self.getVertices()], :]
    
    def getNormal(self):
        return getFaceNormal(self.getVerticesPos())
    
    def getArea(self):
        return getPolygonArea(self.getVerticesPos())
    
    def getCentroid(self):
        V = self.getVerticesPos()
        if V.size == 0:
            return np.array([[0, 0, 0]])
        return np.mean(V, 0)
    
    def getPlane(self):
        return Plane3D(self.mesh.VPos[self.startV.ID, :], self.getNormal())

    #Return the closest point inside this face to P
    def getClosestPoint(self, P):
        #TODO: Finish translating to numpy
        print "TODO: getClosestPoint()"
        #return getClosestPoint([v.pos for v in self.getVertices()], P)            

class MeshEdge(object):
    def __init__(self, mesh, v1, v2, ID):
        self.mesh = mesh
        self.ID = ID
        [self.v1, self.v2] = [v1, v2]
        [self.f1, self.f2] = [None, None]
    
    def vertexAcross(self, startV):
        if startV == self.v1:
            return self.v2
        if startV == self.v2:
            return self.v1
        sys.stderr.write("Warning (vertexAcross): Vertex not member of edge\n")
        return None
    
    def addFace(self, face, v1):
        if self.f1 == None:
            self.f1 = face
        elif self.f2 == None:
            self.f2 = face
        else:
            sys.stderr.write("Cannot add face to edge; already 2 there\n")
    
    #Remove pointer to face
    def removeFace(self, face):
        if self.f1 == face:
            self.f1 = None
        elif self.f2 == face:
            self.f2 = None
        else:
            sys.stderr.write("Cannot remove edge pointer to face that was never part of edge\n")
    
    def faceAcross(self, startF):
        if startF == self.f1:
            return self.f2
        if startF == self.f2:
            return self.f1
        sys.stderr.write("Warning (faceAcross): Face not member of edge\n")
        return None
    
    def getCenter(self):
        V = np.zeros((2, 3))
        V[0, :] = self.mesh.VPos[self.v1.ID, :]
        V[1, :] = self.mesh.VPos[self.v2.ID, :]
        return np.mean(V, 0)
    
    def numAttachedFaces(self):
        ret = 0
        if self.f1:
            ret = ret + 1
        if self.f2:
            ret = ret + 1
        return ret

def getFaceInCommon(e1, e2):
    e2faces = []
    if e2.f1 != None:
        e2faces.append(e2.f1)
    if e2.f2 != None:
        e2faces.append(e2.f2)
    if e1.f1 in e2faces:
        return e1.f1
    if e1.f2 in e2faces:
        return e1.f2
    return None

def getEdgeInCommon(v1, v2):
    for e in v1.edges:
        if e.vertexAcross(v1) is v2:
            return e
    return None

def getVertexInCommon(e1, e2):
    v = [e1.v1, e1.v2, e2.v1, e2.v2]
    for i in range(4):
        for j in range(i+1, 4):
            if v[i] == v[j]:
                return v[i]
    return None

#############################################################
####                    POLYGON MESH                    #####
#############################################################

class PolyMesh(object):
    def __init__(self):
        self.idxpointSize = 20
        self.needsDisplayUpdate = True #Vertex, color, texture, and index buffers need updating
        self.vertices = []
        self.edges = []
        self.faces = []
        #Client side numpy arrays for geometric information
        self.VPos = np.zeros((0, 3)) #Client side vertex buffer
        self.VNormals = np.zeros((0, 3)) #Client side vertex normals
        self.VColors = np.zeros((0, 3)) #Client side color buffer
        self.VTexCoords = np.zeros((0, 2)) #Client side texture coordinate buffer
        self.ITris = np.zeros((0, 3)) #Client side triangle index buffer
        self.EdgeLines = np.zeros((0, 3)) #Client side edge lines
        #Buffer pointers
        self.texID = -1
        self.VPosVBO = None
        self.VNormalsVBO = None
        self.VColorsVBO = None
        self.VTexCoordsVBO = None
        self.IndexVBO = None
        self.EdgeLinesVBO = None
        self.VNormalLinesVBO = None
        self.FNormalLinesVBO = None
        #Display lists
        self.IndexDisplayList = -1
        #Pointers to a representative vertex in different 
        #connected components
        self.components = []
    
    def Clone(self):
        print "TODO"
    
    #Return the edge between v1 and v2 if it exists, or
    #return None if an edge has not yet been created 
    #between them
    def getEdge(self, v1, v2):
        edge = v1.edges & v2.edges
        if len(edge) > 1:
            sys.stderr.write("Warning: More than one edge found on vertex list intersection\n")
        for e in edge:
            return e
        return None
    
    def getVerticesCols(self):
        return self.VPos.T
    
    #############################################################
    ####                ADD/REMOVE METHODS                  #####
    #############################################################

    def addVertex(self, pos, color = np.array([0.5, 0.5, 0.5]), texCoords = np.array([0, 0])):
        vertex = MeshVertex(self, pos, color, texCoords, len(self.vertices))
        self.vertices.append(vertex)
        return vertex
    
    #Create an edge between v1 and v2 and return it
    #This function assumes v1 and v2 are valid vertices in the mesh
    def addEdge(self, v1, v2):
        edge = MeshEdge(self, v1, v2, len(self.edges))
        self.edges.append(edge)
        v1.edges.add(edge)
        v2.edges.add(edge)
        return edge
    
    #Given a list of pointers to mesh vertices in CCW order
    #create a face object from them
    def addFace(self, meshVerts):
        verts = np.array([self.VPos[v.ID, :] for v in meshVerts])
        if not arePlanar(verts):
            sys.stderr.write("Error: Trying to add mesh face that is not planar\n")
            for v in verts:
                print v
            return None
        if not are2DConvex(verts):
            sys.stderr.write("Error: Trying to add mesh face that is not convex\n")
            return None
        face = MeshFace(self, len(self.faces))
        face.startV = meshVerts[0]
        for i in range(0, len(meshVerts)):
            v1 = meshVerts[i]
            v2 = meshVerts[(i+1)%len(meshVerts)]
            edge = self.getEdge(v1, v2)
            if edge == None:
                edge = self.addEdge(v1, v2)
            face.edges.append(edge)
            edge.addFace(face, v1) #Add pointer to face from edge
        self.faces.append(face)
        return face
    
    #Remove the face from the list of faces and remove the pointers
    #from all edges to this face
    def removeFace(self, face):
        #Swap the face to remove with the last face (O(1) removal)
        self.faces[face.ID] = self.faces[-1]
        self.faces[face.ID].ID = face.ID #Update ID of swapped face
        face.ID = -1
        self.faces.pop()
        #Remove pointers from all of the face's edges
        for edge in face.edges:
            edge.removeFace(face)
    
    #Remove this edge from the list of edges and remove 
    #references to the edge from both of its vertices
    #(NOTE: This function is not responsible for cleaning up
    #faces that may have used this edge; that is up to the client)
    def removeEdge(self, edge):
        #Swap the edge to remove with the last edge
        self.edges[edge.ID] = self.edges[-1]
        self.edges[edge.ID].ID = edge.ID #Update ID of swapped face
        edge.ID = -1
        self.edges.pop()
        #Remove pointers from the two vertices that make up this edge
        edge.v1.edges.remove(edge)
        edge.v2.edges.remove(edge)
    
    #Remove this vertex from the list of vertices
    #NOTE: This function is not responsible for cleaning up any of
    #the edges or faces that may have used this vertex
    def removeVertex(self, vertex):
        #Also swap color, vertex, and texture information in buffers
        self.vertices[vertex.ID] = self.vertices[-1]
        self.VPos[vertex.ID, :] = self.VPos[-1, :]
        self.VColors[vertex.ID, :] = self.VColors[-1, :]
        self.VTexCoords[vertex.ID, :] = self.VTexCoords[-1, :]
        self.vertices[vertex.ID].ID = vertex.ID
        vertex.ID = -1
        self.vertices.pop()
    
    #############################################################
    ####    TOPOLOGY SUBDIVISION AND REMESHING METHODS      #####
    #############################################################
    
    #Split every face into K+1 faces by creating a new vertex
    #at the midpoint of every edge, removing the original face, 
    #creating a new face connecting all the new vertices, and
    #creating triangular faces to fill in the gaps
    def splitFaces(self):
        #First subdivide each edge
        for e in self.edges:
            P = e.getCenter()
            e.centerVertex = self.addVertex(P)
        faces = list(self.faces)
        #Now create the new faces
        for f in faces:
            self.removeFace(f)
            #Add the inner face
            fInnerVerts = [e.centerVertex for e in f.edges]
            self.addFace(fInnerVerts)
            #Add the triangles around the border
            innerVerts = [0]*len(f.edges)
            outerVerts = [0]*len(f.edges)
            outerV = f.startV
            for i in range(0, len(f.edges)):
                outerV = f.edges[i].vertexAcross(outerV)
                outerVerts[i] = outerV
                innerVerts[i] = f.edges[i].centerVertex
            for i in range(0, len(f.edges)):
                triVerts = [innerVerts[i], outerVerts[i], innerVerts[(i+1)%len(innerVerts)]]
                self.addFace(triVerts)
        #Remove all edges that were on the original faces
        for f in faces:
            for e in f.edges:
                if e.ID != -1: #If the edge has not already been removed
                    self.removeEdge(e)
        self.needsDisplayUpdate = True
    
    #Split every face into N triangles by creating a vertex at
    #the centroid of each face
    def starRemeshFaces(self, faces):
        for f in faces:
            #TODO: Implement normal meshes this way (move centroid along normal)??
            centroidP = f.getCentroid()
            centroid = self.addVertex(centroidP)
            #TODO: Update texture coordinates properly
            verts = f.getVertices()
            #Remove face and replace with N triangles
            self.removeFace(f)
            for i in range(0, len(verts)):
                v1 = verts[i]
                v2 = verts[(i+1)%len(verts)]
                newVerts = [v1, v2, centroid]
                newFace = self.addFace(newVerts)
        self.needsDisplayUpdate = True
        self.needsIndexDisplayUpdate = True
    
    def starRemesh(self):
        #Copy over the face list since it's about to be modified
        faces = list(self.faces)
        self.starRemeshFaces(faces)
    
    #Triangulate all faces that are not triangular by using
    #the star scheme
    def starTriangulate(self):
        faces = []
        for f in self.faces:
            if len(f.edges) > 3:
                faces.append(f)
        self.starRemeshFaces(faces)
    
    #This function works best starting with triangular meshes
    def evenTriangleRemesh(self):
        for e in self.edges:
            pos = e.getCenter()
            color = 0.5*self.VColors[e.v1.ID, :] + 0.5*self.VColors[e.v2.ID, :]
            texCoords = 0.5*self.VTexCoords[e.v1.ID, :] + 0.5*self.VTexCoords[e.v2.ID, :]
            e.centerVertex = self.addVertex(pos, color, texCoords)
        facesToAdd = []
        for f in self.faces:
            #Add the faces along the outside
            for i in range(0, len(f.edges)):
                e1 = f.edges[i]
                e2 = f.edges[(i+1)%len(f.edges)]
                v = getVertexInCommon(e1, e2)
                facesToAdd.append([e1.centerVertex, v, e2.centerVertex])
            #Add the face in the center
            facesToAdd.append([e.centerVertex for e in f.edges])
        
        #Remove every face and every original edge
        #and add the new faces (which will implicitly 
        #add the new split edges)
        self.faces = []
        self.edges = []
        for v in self.vertices:
            v.edges.clear()
        for f in facesToAdd:
            self.addFace(f)

    #Divide every polygonal face with N sides into (N-2) triangles
    #This assumes the polygonal faces are convex
    def minTrianglesRemesh(self):
        faces = list(self.faces)
        for f in faces:
            verts = f.getVertices()
            if len(verts) <= 3:
                continue
            #Remove face and replace with (N-2) triangles
            self.removeFace(f)
            v0 = verts[0]
            for i in range(1, len(verts)-1):
                v1 = verts[i]
                v2 = verts[i+1]
                newFace = self.addFace([v0, v1, v2])
        self.needsDisplayUpdate = True
    
    def flipNormals(self):
        for f in self.faces:
            f.flipNormal()
    
    def getConnectedComponents(self):
        self.components = []
        for v in self.vertices:
            if v.component == -1:
                self.components.append(v)
                stack = [v]
                while len(stack) > 0:
                    vcurr = stack[-1]
                    if vcurr.component == -1:
                        vcurr.component = len(self.components)-1
                        stack.pop()
                        for vNeighbor in vcurr.getVertexNeighbors():
                            stack.append(vNeighbor)
                    else:
                        stack.pop()
    
    def getConnectedComponentCounts(self):
        counts = [0]*len(self.components)
        for v in self.vertices:
            if v.component > -1:
                counts[v.component] = counts[v.component] + 1
        return counts
    
    def deleteAllButLargestConnectedComponent(self):
        if len(self.components) == 0:
            self.getConnectedComponents()
        counts = self.getConnectedComponentCounts()
        largestComponent = 0
        largestCount = 0
        for i in range(0, len(counts)):
            if counts[i] > largestCount:
                largestCount = counts[i]
                largestComponent = i
        facesToDel = []
        edgesToDel = []
        #Delete faces first, then edges, then vertices
        for f in self.faces:
            if f.startV.component != largestComponent:
                edgesToDel = edgesToDel + f.edges
                facesToDel.append(f)
        for f in facesToDel:
            self.removeFace(f)
        for e in edgesToDel:
            if e.ID != -1:
                self.removeEdge(e)
        verticesToDel = []
        for v in self.vertices:
            if v.component != largestComponent:
                verticesToDel.append(v)
        for v in verticesToDel:
            self.removeVertex(v)
        #Now update the connected components list
        if len(self.vertices) > 0:
            self.components = [self.vertices[0]]
            for v in self.vertices:
                v.component = 0
        self.needsDisplayUpdate = True
    
    #Fill hole with the "advancing front" method
    #but keep it simple for now; not tests for self intersections
    def fillHole(self, hole):
        #TODO: Maintain a min priority queue that indexes based on angle
        #for i in range(0, len(hole)):
        #    vb = hole[(i+len(hole)-1)%len(hole)]
        #    v = hole[i]
        #    va = hole[(i+1)%len(hole)]
        #    v.before = vb
        #    v.after = va
        c = [0, 0, 0]
        for mV in hole:
            c = c + self.VPos[mV.ID, :]
        c = (1.0/len(hole))*c
        vCenter = self.addVertex(c)
        for i in range(0, len(hole)):
            v1 = hole[i]
            v2 = hole[(i+1)%len(hole)]
            self.addFace([vCenter, v1, v2])
        
    
    def fillHoles(self, slicedHolesOnly = False):
        holes = []
        origEdges = self.edges[:]
        for e in origEdges:
            if e.numAttachedFaces() == 1 and ((not slicedHolesOnly) or e.v1.borderVertex):
                loop = [e.v1, e.v2]
                finished = False
                while not finished:
                    foundNext = False
                    for v in loop[-1].getVertexNeighbors():
                        if v is loop[-2]:
                            #Make sure it doesn't accidentally back up
                            continue
                        elif v is loop[0]:
                            #It's looped back on itself so we're done
                            finished = True
                        else:
                            e = getEdgeInCommon(loop[-1], v)
                            if not e:
                                print "Warning: Edge not found in common while trying to trace hole boundary"
                                finished = True
                                break
                            elif e.numAttachedFaces() == 1:
                                foundNext = True
                                loop.append(v)
                                break
                    if not foundNext and not finished:
                        print "Warning: Unable to close hole"
                        break
                print "Found hole of size %i"%len(loop)
                self.fillHole(loop)
        self.needsDisplayUpdate = True
        self.needsIndexDisplayUpdate = True

    #############################################################
    ####                 GEOMETRY METHODS                   #####
    #############################################################

    #Transformations are simple because geometry information is only
    #stored in the vertices
    def Transform(self, matrix):
        self.VPos = matrix.dot(self.VPos.T).T
    
    def Translate(self, dV):
        self.VPos = self.VPos + np.reshape(dV, [1, 3])
    
    def Scale(self, dx, dy, dz):
        S = np.zeros((1, 3))
        S[0, :] = [dx, dy, dz]
        self.VPos = S*self.VPos

    def getCentroid(self):
        return np.mean(self.VPos, 0)
    
    def getBBox(self):
        if self.VPos.shape[0] == 0:
            print "Warning: PolyMesh.getBBox(): Adding bbox but no vertices"
            return BBox3D()
        bbox = BBox3D()
        bbox.fromPoints(self.VPos)
        return bbox
    
    #Use PCA to find the principal axes of the vertices
    def getPrincipalAxes(self):
        X = self.VPos - self.getCentroid()
        XTX = (X.T).dot(X)
        (lambdas, axes) = linalg.eig(XTX)
        #Put the eigenvalues in decreasing order
        idx = lambdas.argsort()[::-1]
        lambdas = lambdas[idx]
        axes = axes[:, idx]
        T = X.dot(axes)
        maxProj = T.max(0)
        minProj = T.min(0)
        axes = axes.T #Put each axis on each row to be consistent with everything else
        return (axes, maxProj, minProj)        
    
    #Delete the parts of the mesh below "plane".  If fillHoles
    #is true, plug up the holes that result from the cut
    def sliceBelowPlane(self, plane, fillHoles = True):
        facesToDel = []
        edgesToDel = []
        verticesToDel = []
        facesToAdd = []
        newVertices = []
        borderVertex = None
        for e in self.edges:
            #Keep track of edges which intersect the plane
            #and the vertex that represents that intersection
            e.centerVertex = None
        for v in self.vertices:
            #Keep track of which vertices are introduced at the plane slice
            v.borderVertex = False
        for f in self.faces:
            v1 = f.startV
            deleteFace = False
            splitFaceStartE = None
            splitFaceStartEIndex = -1
            splitFaceStartV = None
            splitFaceEndE = None
            for i in range(0, len(f.edges)):
                e = f.edges[i]
                v2 = e.vertexAcross(v1)
                distV1 = plane.distFromPlane(self.VPos[v1.ID, :])
                distV2 = plane.distFromPlane(self.VPos[v2.ID, :])
                #Mark all vertices below the plane for deletion
                #(this will count vertices more than once but 
                #because mesh is manifold the amortized size will
                #still be linear in the number of vertices)
                if distV1 < 0:
                    verticesToDel.append(v1)
                #Mark all edges that are fully or partially below
                #the plane for deletion.  This list will hold
                #every such edge at most twice since 2 faces meet
                #at an edge in a manifold mesh
                if distV1 < 0 or distV2 < 0:
                    edgesToDel.append(e)
                    deleteFace = True
                #Edge goes from negative side to plus side of plane
                if distV1 < 0 and distV2 >= 0:
                    if not e.centerVertex:
                        line = Line3D(self.VPos[v1.ID, :], self.VPos[v2.ID, :] - self.VPos[v1.ID, :])
                        [t, P] = line.intersectPlane(plane)
                        if P:
                            newColor = (1-t)*self.VColors[v1.ID, :] + t*self.VColors[v2.ID, :]
                            newTexCoords = (1-t)*self.VTexCoords[v1.ID, :] + t*self.VTexCoords[v2.ID, :]
                            e.centerVertex = self.addVertex(P, newColor, newTexCoord)
                            e.centerVertex.borderVertex = True
                            borderVertex = e.centerVertex
                    if e.centerVertex:
                        splitFaceStartEIndex = i
                        splitFaceStartV = v1
                        splitFaceStartE = e
                #Edge goes from plus side to negative side of plane
                if distV1 >= 0 and distV2 <= 0:
                    if not e.centerVertex:
                        line = Line3D(self.VPos[v1.ID, :], self.VPos[v2.ID, :] - self.VPos[v1.ID, :])
                        [t, P] = line.intersectPlane(plane)
                        if P:
                            newColor = (1-t)*self.VColors[v1.ID, :] + t*self.VColors[v2.ID, :]
                            newTexCoords = (1-t)*self.VTexCoords[v1.ID, :] + t*self.VTexCoords[v2.ID, :]
                            e.centerVertex = self.addVertex(P, newColor, newTexCoords)
                            e.centerVertex.texCoords = newTexCoords
                            e.centerVertex.borderVertex = True
                            borderVertex = e.centerVertex
                    if e.centerVertex:
                        splitFaceEndE = e                    
                v1 = v2
            if deleteFace:
                facesToDel.append(f)
            #Walk along the split part of the face on the positive
            #side of the plane
            if splitFaceStartE and splitFaceEndE:
                newFace = [splitFaceStartE.centerVertex]
                i = splitFaceStartEIndex
                e = splitFaceStartE
                v1 = splitFaceStartV
                while e != splitFaceEndE:
                    v1 = e.vertexAcross(v1)
                    newFace.append(v1)
                    i = (i+1)%len(f.edges)
                    e = f.edges[i]
                newFace.append(splitFaceEndE.centerVertex)
                facesToAdd.append(newFace)
        #First remove all faces that are no longer relevant
        for f in facesToDel:
            self.removeFace(f)    
        #Now remove edges that are no longer relevant
        for e in edgesToDel:
            if e.ID != -1:
                self.removeEdge(e)
        #Now remove vertices that are no longer relevant
        for v in verticesToDel:
            if v.ID != -1:
                self.removeVertex(v)
        #Add new faces
        for f in facesToAdd:
            self.addFace(f)
        if fillHoles:
            self.fillHoles(slicedHolesOnly = True)
        self.needsDisplayUpdate = True
        self.needsIndexDisplayUpdate = True
    
    def sliceAbovePlane(self, plane, fillHoles = True):
        planeNeg = Plane3D(plane.P0, plane.N)
        planeNeg.initFromEquation(-plane.A, -plane.B, -plane.C, -plane.D)
        self.sliceBelowPlane(planeNeg, fillHoles)
    
    def flipAcrossPlane(self, plane):
        P0 = plane.P0
        N = plane.N
        for V in self.vertices:
            P = self.VPos[V.ID, :]
            dP = P - P0
            dPPar = projVec(dP, N)
            dPPerp = dP - dPPar
            self.VPos[V.ID, :] = P0 - dPPar + dPPerp
        self.needsDisplayUpdate = True
    
    #Uniformly sample points on this mesh at random, taking into consideration
    #the area of triangle faces
    #Return the points and normal estimates
    #colPoints: Whether to return the points/normals in columns or rows
    def randomlySamplePoints(self, NPoints, colPoints = True):
        if self.needsDisplayUpdate:
            #Make sure the triangle buffer is in place even if the rest of
            #the buffers haven't been setup yet
            self.updateTris()
        ###Step 1: Compute cross product of all face triangles and use to compute
        #areas and normals (very similar to code used to compute vertex normals)
        
        #Vectors spanning two triangle edges
        P0 = self.VPos[self.ITris[:, 0], :]
        P1 = self.VPos[self.ITris[:, 1], :]
        P2 = self.VPos[self.ITris[:, 2], :]
        V1 = P1 - P0
        V2 = P2 - P0
        FNormals = np.cross(V1, V2)
        FAreas = np.sqrt(np.sum(FNormals**2, 1)).flatten()
        
        #Get rid of zero area faces and update points
        self.ITris = self.ITris[FAreas > 0, :]
        FNormals = FNormals[FAreas > 0, :]
        FAreas = FAreas[FAreas > 0]
        P0 = self.VPos[self.ITris[:, 0], :]
        P1 = self.VPos[self.ITris[:, 1], :]
        P2 = self.VPos[self.ITris[:, 2], :]
        
        #Compute normals
        NTris = self.ITris.shape[0]
        FNormals = FNormals/FAreas[:, None]
        FAreas = 0.5*FAreas
        self.FNormals = FNormals
        self.FCentroid = 0*FNormals
        self.VNormals = 0*self.VPos
        VAreas = np.zeros(self.VPos.shape[0])
        for k in range(3):
            self.VNormals[self.ITris[:, k], :] += FAreas[:, None]*FNormals
            VAreas[self.ITris[:, k]] += FAreas
        #Normalize normals
        VAreas[VAreas == 0] = 1
        self.VNormals = self.VNormals / VAreas[:, None]
        
        ###Step 2: Randomly sample points based on areas
        FAreas = FAreas/np.sum(FAreas)
        AreasC = np.cumsum(FAreas)
        samples = np.sort(np.random.rand(NPoints))
        #Figure out how many samples there are for each face
        FSamples = np.zeros(NTris)
        fidx = 0
        for s in samples:
            while s > AreasC[fidx]:
                fidx += 1
            FSamples[fidx] += 1
        #Now initialize an array that stores the triangle sample indices
        tidx = np.zeros(NPoints, dtype=np.int64)
        idx = 0
        for i in range(len(FSamples)):
            tidx[idx:idx+FSamples[i]] = i
            idx += FSamples[i]
        N = np.zeros((NPoints, 3)) #Allocate space for normals
        idx = 0
        
        #Vector used to determine if points need to be flipped across parallelogram
        V3 = P2 - P1
        V3 = V3/np.sqrt(np.sum(V3**2, 1))[:, None] #Normalize
        
        #Randomly sample points on each face        
        #Generate random points uniformly in parallelogram
        u = np.random.rand(NPoints, 1)
        v = np.random.rand(NPoints, 1)
        Ps = u*V1[tidx, :] + P0[tidx, :]
        Ps += v*V2[tidx, :]
        #Flip over points which are on the other side of the triangle
        dP = Ps - P1[tidx, :]
        proj = np.sum(dP*V3[tidx, :], 1)
        dPPar = V3[tidx, :]*proj[:, None] #Parallel project onto edge
        dPPerp = dP - dPPar
        Qs = Ps - dPPerp
        dP0QSqr = np.sum((Qs - P0[tidx, :])**2, 1)
        dP0PSqr = np.sum((Ps - P0[tidx, :])**2, 1)
        idxreg = np.arange(NPoints, dtype=np.int64)
        idxflip = idxreg[dP0QSqr < dP0PSqr]
        u[idxflip, :] = 1 - u[idxflip, :]
        v[idxflip, :] = 1 - v[idxflip, :]
        Ps[idxflip, :] = P0[tidx[idxflip], :] + u[idxflip, :]*V1[tidx[idxflip], :] + v[idxflip, :]*V2[tidx[idxflip], :]
        
        #Step 3: Compute normals of sampled points by barycentric interpolation
        Ns = u*self.VNormals[self.ITris[tidx, 1], :]
        Ns += v*self.VNormals[self.ITris[tidx, 2], :] 
        Ns += (1-u-v)*self.VNormals[self.ITris[tidx, 0], :]
        
        if colPoints:
            return (Ps.T, Ns.T)
        return (Ps, Ns)
    
    #############################################################
    ####                INPUT/OUTPUT METHODS                #####
    #############################################################
    def loadFile(self, filename):
        suffix = re.split("\.", filename)[-1]
        if suffix == "off":
            self.loadOffFile(filename)
        elif suffix == "toff":
            self.loadTOffFile(filename)
        elif suffix == "obj":
            self.loadObjFile(filename)
        else:
            print "Unsupported file suffix (%s) for loading mesh"%(suffix, filename)
        self.needsDisplayUpdate = True
        self.needsIndexDisplayUpdate = True
    
    def saveFile(self, filename, verbose = False):
        suffix = re.split("\.", filename)[-1]
        if suffix == "off":
            self.saveOffFile(filename, verbose)
        elif suffix == "obj":
            self.saveObjFile(filename, verbose)
        elif suffix == "ply":
            self.savePlyFile(filename, verbose)
        else:
            print "Unsupported file suffix (%s) for saving mesh %s"%(suffix, filename)        
    
    def loadOffFile(self, filename):
        fin = open(filename, 'r')
        nVertices = 0
        nFaces = 0
        nEdges = 0
        lineCount = 0
        face = 0
        vertex = 0
        divideColor = False
        for line in fin:
            lineCount = lineCount+1
            fields = line.split() #Splits whitespace by default
            if len(fields) == 0: #Blank line
                continue
            if fields[0][0] in ['#', '\0', ' '] or len(fields[0]) == 0:
                continue
            #Check section
            if nVertices == 0:
                if fields[0] == "OFF" or fields[0] == "COFF":
                    if len(fields) > 2:
                        fields[1:4] = [int(field) for field in fields]
                        [nVertices, nFaces, nEdges] = fields[1:4]  
                        #Pre-allocate vertex arrays    
                        self.VPos = np.zeros((nVertices, 3)) 
                        self.VColors = np.zeros((nVertices, 3))
                        self.VTexCoords = np.zeros((nVertices, 2))
                    if fields[0] == "COFF":
                        divideColor = True            
                else:
                    fields[0:3] = [int(field) for field in fields]
                    [nVertices, nFaces, nEdges] = fields[0:3]
            elif vertex < nVertices:
                fields = [float(i) for i in fields]
                P = [fields[0],fields[1], fields[2]]
                color = np.array([0.5, 0.5, 0.5]) #Gray by default
                if len(fields) >= 6:
                    #There is color information
                    if divideColor:
                        color = [float(c)/255.0 for c in fields[3:6]]
                    else:
                        color = [float(c) for c in fields[3:6]]
                self.addVertex(P, color)
                vertex = vertex+1
            elif face < nFaces:
                #Assume the vertices are specified in CCW order
                fields = [int(i) for i in fields]
                meshVerts = fields[1:fields[0]+1]
                verts = [self.vertices[i] for i in meshVerts]
                self.addFace(verts)
                face = face+1
        fin.close()
        if np.max(self.VColors) > 1:
            #Rescale colors
            self.VColors = self.VColors / 255.0
    
    #My own "TOFF" format, which is like OFF with texture
    def loadTOffFile(self, filename):
        fin = open(filename, 'r')
        nVertices = 0
        nFaces = 0
        nEdges = 0
        lineCount = 0
        face = 0
        vertex = 0
        textureName = ""
        for line in fin:
            lineCount = lineCount+1
            fields = line.split() #Splits whitespace by default
            if len(fields) == 0: #Blank line
                continue
            if fields[0][0] in ['#', '\0', ' '] or len(fields[0]) == 0:
                continue
            #Check section
            if nVertices == 0:
                if fields[0] == "TOFF":
                    textureName = fields[1]
                    self.texID = loadTexture(textureName)    
                else:
                    fields[0:3] = [int(field) for field in fields]
                    [nVertices, nFaces, nEdges] = fields[0:3]
                    #Pre-allocate vertex arrays    
                    self.VPos = np.zeros((nVertices, 3)) 
                    self.VColors = np.zeros((nVertices, 3))
                    self.VTexCoords = np.zeros((nVertices, 2))
            elif vertex < nVertices:
                fields = [float(i) for i in fields]
                P = [fields[0],fields[1], fields[2]]
                v = self.addVertex(P)
                v.texCoords = [fields[3], fields[4]]
                vertex = vertex+1
            elif face < nFaces:
                #Assume the vertices are specified in CCW order
                fields = [int(i) for i in fields]
                meshVerts = fields[1:fields[0]+1]
                verts = [self.vertices[i] for i in meshVerts]
                self.addFace(verts)
                face = face+1
        fin.close()
            
    def saveOffFile(self, filename, verbose = False, outputColors = True, output255 = False):
        nV = len(self.vertices)
        nE = len(self.edges)
        nF = len(self.faces)
        fout = open(filename, "w")
        #fout.write("#Generated with Chris Tralie's G-RFLCT Library\n")
        #fout.write("#http://www.github.com/ctralie/G-RFLCT\n")
        fout.write("COFF\n%i %i %i\n"%(nV, nF, 0))
        for v in self.vertices:
            fout.write("%g %g %g"%tuple(self.VPos[v.ID, :]))
            if outputColors:
                c = self.VColors[v.ID, :]
                if output255:
                    fout.write(" %i %i %i"%tuple(np.round(c)))
                else:
                    fout.write(" %g %g %g"%tuple(c))
            fout.write("\n")
        for f in self.faces:
            verts = f.getVertices()
            fout.write("%i "%(len(verts)))
            for v in verts:
                fout.write("%i "%(v.ID))
            fout.write("\n")
        fout.close()
        if verbose:
            print "Saved file to %s"%filename

    def savePlyFile(self, filename, verbose = False, outputColors = True, output255 = True):
        nV = len(self.vertices)
        nE = len(self.edges)
        nF = len(self.faces)
        fout = open(filename, "w")
        fout.write("ply\nformat ascii 1.0\nelement vertex %i\n"%nV)
        fout.write("property float x\nproperty float y\nproperty float z\n")
        if outputColors:
            fout.write("property uchar red\nproperty uchar green\nproperty uchar blue\n")
        fout.write("element face %i\n"%nF)
        fout.write("property list uchar int vertex_indices\nend_header\n")
        for v in self.vertices:
            fout.write("%g %g %g"%tuple(self.VPos[v.ID, :]))
            if outputColors:
                c = self.VColors[v.ID, :]
                if output255:
                    fout.write(" %i %i %i"%tuple(np.round(c)))
                else:
                    fout.write(" %g %g %g"%tuple(c))
            fout.write("\n")
        for f in self.faces:
            verts = f.getVertices()
            fout.write("%i "%(len(verts)))
            for v in verts:
                fout.write("%i "%(v.ID))
            fout.write("\n")
        fout.close()
        if verbose:
            print "Saved file to %s"%filename
        
    def loadObjFile(self, filename):
        #TODO: Right now vertex normals, face normals, and texture coordinates are ignored
        #Later incorporate them??
        fin = open(filename, 'r')
        for line in fin:
            fields = line.split()
            if len(fields) == 0: #Blank line
                continue
            if fields[0][0] in ['#', '\0', ' '] or len(fields[0]) == 0:
                continue
            if fields[0] == "v":
                coords = [float(i) for i in fields[1:4]]
                self.addVertex([coords[0], coords[1], coords[2]])
            if fields[0] == "f":
                #Indices are numbered starting at 1 (so need to subtract that off)
                indices = [int(re.split("/",s)[0])-1 for s in fields[1:]]
                verts = [self.vertices[i] for i in indices]
                self.addFace(verts)
        fin.close()
    
    def saveObjFile(self, filename, verbose = False):
        fout = open(filename, "w")
        fout.write("#Generated with Chris Tralie's G-RFLCT Library\n")
        fout.write("#http://www.github.com/ctralie/G-RFLCT\n")
        fout.flush()
        for v in self.vertices:
            fout.write("%g %g %g"%tuple(self.VPos[v.ID, :]))
        for f in self.faces:
            verts = f.getVertices()
            fout.write("f ")
            i = 0
            for i in range(0, len(verts)):
                v = verts[i]
                fout.write("%i"%(v.ID+1))#Indices are numbered starting at 1
                if i < len(verts) - 1:
                    fout.write(" ")
            fout.write("\n")
        fout.close()
        if verbose:
            print "Saved file to %s"%filename

    #############################################################
    ####                     RENDERING                      #####
    #############################################################
    
    #Figure out triangles that go into the index buffer, splitting faces
    #with more than three vertices into triangles (assuming convexity)
    def updateTris(self):
        NTris = np.sum(np.array([len(f.edges)-2 for f in self.faces]))
        self.ITris = np.zeros((NTris, 3), dtype = np.int32)
        idx = 0
        for f in self.faces:
            IDs = [v.ID for v in f.getVertices()]
            for k in range(len(IDs)-2):
                self.ITris[idx, :] = [IDs[0]] + IDs[k+1:k+3]
                idx += 1
    
    def updateNormalBuffer(self):
        #First compute cross products of all face triangles
        V1 = self.VPos[self.ITris[:, 1], :] - self.VPos[self.ITris[:, 0], :]
        V2 = self.VPos[self.ITris[:, 2], :] - self.VPos[self.ITris[:, 0], :]
        FNormals = np.cross(V1, V2)
        FAreas = np.reshape(np.sqrt(np.sum(FNormals**2, 1)), (FNormals.shape[0], 1))
        FAreas[FAreas == 0] = 1
        FNormals = FNormals/FAreas
        self.FNormals = FNormals
        self.FCentroid = 0*FNormals
        self.VNormals = 0*self.VPos
        VAreas = np.zeros((self.VPos.shape[0], 1))
        for k in range(3):
            self.VNormals[self.ITris[:, k], :] += FAreas*FNormals
            VAreas[self.ITris[:, k]] += FAreas
            self.FCentroid += self.VPos[self.ITris[:, k], :]
        self.FCentroid /= 3
        #Normalize normals
        VAreas[VAreas == 0] = 1
        self.VNormals = self.VNormals / VAreas
    
    #buffersOnly: True if there is no mesh structure other than VPos
    #VColors and ITris
    def performDisplayUpdate(self, buffersOnly = False):
        #Clear the buffers for the invalidated mesh
        if self.VPosVBO:
            self.VPosVBO.delete()
        if self.VNormalsVBO:
            self.VNormalsVBO.delete()
        if self.VColorsVBO:
            self.VColorsVBO.delete()
        if self.VTexCoordsVBO:
            self.VTexCoordsVBO.delete()
        if self.IndexVBO:
            self.IndexVBO.delete()
        if self.EdgeLinesVBO:
            self.EdgeLinesVBO.delete()
        if self.VNormalLinesVBO:
            self.VNormalLinesVBO.delete()
        if self.FNormalLinesVBO:
            self.FNormalLinesVBO.delete()
        self.VPosVBO = vbo.VBO(np.array(self.VPos, dtype=np.float32))
        self.VColorsVBO = vbo.VBO(np.array(self.VColors, dtype=np.float32))
        self.VTexCoordsVBO = vbo.VBO(np.array(self.VTexCoords, dtype=np.float32))
        if not buffersOnly:
            self.updateTris()
        self.IndexVBO = vbo.VBO(self.ITris, target=GL_ELEMENT_ARRAY_BUFFER)
        
        #Update edges buffer
        if buffersOnly:
            #Use triangle faces to add edges (will be redundancy but is faster
            #tradeoff for rendering only)
            NTris = self.ITris.shape[0]
            self.EdgeLines = np.zeros((NTris*3*2, 3))
            for k in range(3):
                istart = k*NTris*2
                self.EdgeLines[istart:istart+NTris*2:2, :] = self.VPos[self.ITris[:, k], :]
                self.EdgeLines[istart+1:istart+NTris*2:2, :] = self.VPos[self.ITris[:, (k+1)%3], :]
        else:
            #Use the mesh structure to find the edges
            self.EdgeLines = np.zeros((len(self.edges)*2, 3))
            for i in range(len(self.edges)):
                self.EdgeLines[i*2, :] = self.VPos[self.edges[i].v1.ID, :]
                self.EdgeLines[i*2+1, :] = self.VPos[self.edges[i].v2.ID, :]
        self.EdgeLinesVBO = vbo.VBO(np.array(self.EdgeLines, dtype=np.float32))
        
        #Update face and vertex normals
        scale = 0.01*self.getBBox().getDiagLength()
        self.updateNormalBuffer()
        self.VNormalsVBO = vbo.VBO(np.array(self.VNormals, dtype=np.float32))
        VNList = np.zeros((self.VPos.shape[0]*2, 3))
        VNList[np.arange(0, VNList.shape[0], 2), :] = self.VPos
        VNList[np.arange(1, VNList.shape[0], 2), :] = self.VPos + scale*self.VNormals
        self.VNormalLinesVBO = vbo.VBO(np.array(VNList, dtype=np.float32))
        VFList = np.zeros((self.ITris.shape[0]*2, 3))
        VFList[np.arange(0, VFList.shape[0], 2), :] = self.FCentroid
        VFList[np.arange(1, VFList.shape[0], 2), :] = self.FCentroid + scale*self.FNormals
        self.FNormalLinesVBO = vbo.VBO(np.array(VFList, dtype=np.float32))
        
        self.needsDisplayUpdate = False
    
    #vertexColors is an Nx3 numpy array, where N is the number of vertices
    def renderGL(self, drawEdges = False, drawVerts = False, drawFaces = True, drawVertexNormals = True, drawFaceNormals = True, useLighting = True, useTexture = False):
        if self.needsDisplayUpdate:
            self.performDisplayUpdate()
            self.needsDisplayUpdate = False
        #Draw the requested objects
        if useLighting:
            glEnable(GL_LIGHTING)
        else:
            glDisable(GL_LIGHTING)
        if drawFaces:
            glEnableClientState(GL_VERTEX_ARRAY)
            glEnableClientState(GL_COLOR_ARRAY)
            glEnableClientState(GL_NORMAL_ARRAY)
            self.VPosVBO.bind()
            glVertexPointerf(self.VPosVBO)
            self.VNormalsVBO.bind()
            glNormalPointerf(self.VNormalsVBO)
            self.VColorsVBO.bind()
            glColorPointerf(self.VColorsVBO)
            if useTexture and self.VTexCoordsVBO:
                glEnable(GL_TEXTURE_2D)
                glEnableClientState(GL_TEXTURE_COORD_ARRAY)
                self.VTexCoordsVBO.bind()
                glTexCoordPointerf(self.VTexCoordsVBO)
                glBindTexture(GL_TEXTURE_2D, self.texID)
            else:
                glEnable(GL_COLOR_MATERIAL)
                glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE)
            self.IndexVBO.bind()
            glDrawElements(GL_TRIANGLES, 3*self.ITris.shape[0], GL_UNSIGNED_INT, None)
            self.IndexVBO.unbind()
            if useTexture and self.VTexCoordsVBO:
                self.VTexCoordsVBO.unbind()
                glDisableClientState(GL_TEXTURE_COORD_ARRAY)
            self.VPosVBO.unbind()
            self.VNormalsVBO.unbind()
            self.VColorsVBO.unbind()
            glDisableClientState(GL_NORMAL_ARRAY)
            glDisableClientState(GL_COLOR_ARRAY)
            glDisableClientState(GL_VERTEX_ARRAY)
        
        if drawVerts:
            glEnableClientState(GL_VERTEX_ARRAY)
            self.VPosVBO.bind()
            glVertexPointerf(self.VPosVBO)
            glDisable(GL_LIGHTING)
            glPointSize(POINT_SIZE)
            glColor3f(1.0, 0, 0)
            glDrawArrays(GL_POINTS, 0, self.VPos.shape[0])
            self.VPosVBO.unbind()
            glDisableClientState(GL_VERTEX_ARRAY)
        
        if drawEdges:
            glEnableClientState(GL_VERTEX_ARRAY)
            self.EdgeLinesVBO.bind()
            glVertexPointerf(self.EdgeLinesVBO)
            glDisable(GL_LIGHTING)
            glLineWidth(2)
            glColor3f(1.0, 1.0, 0.0)
            glDrawArrays(GL_LINES, 0, self.EdgeLines.shape[0])
            self.EdgeLinesVBO.unbind()
            glDisableClientState(GL_VERTEX_ARRAY)
        
        if drawVertexNormals:
            glEnableClientState(GL_VERTEX_ARRAY)
            self.VNormalLinesVBO.bind()
            glVertexPointerf(self.VNormalLinesVBO)
            glDisable(GL_LIGHTING)
            glColor3f(1.0, 0, 1.0)
            glDrawArrays(GL_LINES, 0, self.VPos.shape[0]*2)
            self.VNormalLinesVBO.unbind()
            glDisableClientState(GL_VERTEX_ARRAY)
        
        if drawFaceNormals:
            glEnableClientState(GL_VERTEX_ARRAY)
            self.FNormalLinesVBO.bind()
            glVertexPointerf(self.VNormalLinesVBO)
            glDisable(GL_LIGHTING)
            glColor3f(1.0, 0, 0)
            glDrawArrays(GL_LINES, 0, self.ITris.shape[0]*2)
            self.FNormalLinesVBO.unbind()
            glDisableClientState(GL_VERTEX_ARRAY)
        
        if useLighting:
            glEnable(GL_LIGHTING)

    #Render the vertices of the mesh as points with colors equal to their indices
    #in the vertex list.  Used to help with vertex selection (will work as long as
    #there are fewer than 2^32 vertices)
    #TODO: Make this faster with vertex buffer instead of display list
    #(but that would require some annoying numpy bitwhacking)
    def renderGLIndices(self):
        if self.needsDisplayUpdate:
            self.performDisplayUpdate()
            self.needsDisplayUpdate = False
        glCallList(self.IndexDisplayList)
        
    def updateIndexDisplayList(self):
        if self.IndexDisplayList != -1: #Deallocate previous display list
            glDeleteLists(self.IndexDisplayList, 1)
        self.IndexDisplayList = glGenLists(1)
        print "Updating index display list"
        glNewList(self.IndexDisplayList, GL_COMPILE)
        glDisable(GL_LIGHTING)
        N = len(self.vertices)
        #First draw all of the faces with index N+1 so that occlusion is
        #taken into proper consideration
        [R, G, B, A] = splitIntoRGBA(N+2)
        glColor4ub(R, G, B, A)
        glEnableClientState(GL_VERTEX_ARRAY)
        self.VPosVBO.bind()
        glVertexPointerf(self.VPosVBO)
        self.IndexVBO.bind()
        glDrawElements(GL_TRIANGLES, 3*self.ITris.shape[0], GL_UNSIGNED_INT, None)
        self.IndexVBO.unbind()
        self.VPosVBO.unbind()
        glPointSize(POINT_SIZE)
        glBegin(GL_POINTS)
        for i in range(0, N):
            P = self.VPos[i, :]
            [R, G, B, A] = splitIntoRGBA(i+1)
            glColor4ub(R, G, B, A)
            glVertex3f(P[0], P[1], P[2])
        glEnd()
        glEnable(GL_LIGHTING)
        glEndList()

    #Slow version with no spatial subdivision
    def getRayIntersection(self, ray):
        t = float("inf")
        Point = None
        Face = None
        for f in self.faces:
            intersection = ray.intersectMeshFace(f)
            if intersection != None:
                if intersection[0] < t:
                    t = intersection[0]
                    Point = intersection[1]
                    Face = f
        return [t, Point, Face]
        return None
    
    def __str__(self):
        nV = len(self.vertices)
        nE = len(self.edges)
        nF = len(self.faces)
        euler = nV-nE+nF
        return "PolyMesh Object: NVertices = %i, NEdges = %i, NFaces = %i, euler=%i"%(nV, nE, nF, euler)    


#############################################################
####               UTILITY FUNCTIONS                    #####
#############################################################
def saveOffFileExternal(filename, VPos, VColors, ITris):
    #Save off file given buffers, not necessarily in the PolyMesh object
    nV = VPos.shape[0]
    nF = ITris.shape[0]
    fout = open(filename, "w")
    if VColors.size == 0:
        fout.write("OFF\n%i %i %i\n"%(nV, nF, 0))
    else:
        fout.write("COFF\n%i %i %i\n"%(nV, nF, 0))
    for i in range(nV):
        fout.write("%g %g %g"%tuple(VPos[i, :]))
        if VColors.size > 0:
            fout.write(" %g %g %g"%tuple(VColors[i, :]))
        fout.write("\n")
    for i in range(nF):
        fout.write("3 %i %i %i\n"%tuple(ITris[i, :]))
    fout.close()

#Return VPos, VColors, and ITris without creating any structure
#(Assumes triangle mesh)
def loadOffFileExternal(filename):
    fin = open(filename, 'r')
    nVertices = 0
    nFaces = 0
    lineCount = 0
    face = 0
    vertex = 0
    divideColor = False
    VPos = np.zeros((0, 3))
    VColors = np.zeros((0, 3))
    ITris = np.zeros((0, 3))
    for line in fin:
        lineCount = lineCount+1
        fields = line.split() #Splits whitespace by default
        if len(fields) == 0: #Blank line
            continue
        if fields[0][0] in ['#', '\0', ' '] or len(fields[0]) == 0:
            continue
        #Check section
        if nVertices == 0:
            if fields[0] == "OFF" or fields[0] == "COFF":
                if len(fields) > 2:
                    fields[1:4] = [int(field) for field in fields]
                    [nVertices, nFaces, nEdges] = fields[1:4]  
                    print "nVertices = %i, nFaces = %i"%(nVertices, nFaces)
                    #Pre-allocate vertex arrays    
                    VPos = np.zeros((nVertices, 3)) 
                    VColors = np.zeros((nVertices, 3))
                    ITris = np.zeros((nFaces, 3))
                if fields[0] == "COFF":
                    divideColor = True            
            else:
                fields[0:3] = [int(field) for field in fields]
                [nVertices, nFaces, nEdges] = fields[0:3]
                VPos = np.zeros((nVertices, 3)) 
                VColors = np.zeros((nVertices, 3))
                ITris = np.zeros((nFaces, 3))
        elif vertex < nVertices:
            fields = [float(i) for i in fields]
            P = [fields[0],fields[1], fields[2]]
            color = np.array([0.5, 0.5, 0.5]) #Gray by default
            if len(fields) >= 6:
                #There is color information
                if divideColor:
                    color = [float(c)/255.0 for c in fields[3:6]]
                else:
                    color = [float(c) for c in fields[3:6]]
            VPos[vertex, :] = P
            VColors[vertex, :] = color
            vertex = vertex+1
        elif face < nFaces:
            #Assume the vertices are specified in CCW order
            fields = [int(i) for i in fields]
            ITris[face, :] = fields[1:fields[0]+1]
            face = face+1
    fin.close()
    VPos = np.array(VPos, np.float64)
    VColors = np.array(VColors, np.float64)
    ITris = np.array(ITris, np.int32)
    return (VPos, VColors, ITris) 

#Make a surface of revolution around the y-axis (x = 0, z = 0)
#X: An ordrered N x 2 list of points in the XY plane that make up a curve
#NSteps: Number of points to sample around the circle of revolution
def makeSurfaceOfRevolution(X, NSteps):
    N = X.shape[0]
    thetas = np.linspace(0, 2*np.pi, NSteps+1)[0:NSteps]
    M = N*NSteps #Total number of vertices
    VPos = np.zeros((M, 3))
    #First make all of the points
    for i in range(NSteps):
        VPos[N*i:N*(i+1), 1] = X[:, 1] #The Y positions stay the same; it's only XZ that change
        VPos[N*i:N*(i+1), 0] = X[:, 0]*np.cos(thetas[i])
        VPos[N*i:N*(i+1), 2] = X[:, 0]*np.sin(thetas[i])
    #Now make all of the triangle faces connecting them
    L = (N-1)*2 #Number of triangles added for each curve position
    ITris = np.zeros((L*NSteps, 3))
    for i in range(NSteps):
        j = (i+1)%NSteps
        #Upper triangles, CCW order
        idx = np.arange(L*i, L*(i+1), 2)
        ITris[idx, 0] = N*i + np.arange(N-1)
        ITris[idx, 1] = N*j + np.arange(N-1)
        ITris[idx, 2] = N*j + 1 + np.arange(N-1)
        #Lower triangles, CCW order
        idx = np.arange(L*i+1, L*(i+1), 2)
        ITris[idx, 0] = N*i + np.arange(N-1)
        ITris[idx, 1] = N*j + 1 + np.arange(N-1)
        ITris[idx, 2] = N*i + 1 + np.arange(N-1)
    return (VPos, ITris)

#############################################################
####               STANDARD MESHES                      #####
#############################################################

#Helper function for getBoxMesh and addFaceTiles
def makeBoxEdge(mesh, v1, v2, stepSize):
    if stepSize < 0:
        return [v1, v2]
    verts = [v1]
    direc = mesh.VPos[v2.ID, :] - mesh.VPos[v1.ID, :]
    frac = stepSize/np.sqrt(direc.dot(direc))
    #Round to the nearest integer number of tiles
    N = int(math.floor(1.0/frac+0.5))
    if N == 0:
        N = 1
    frac = 1.0/float(N)
    for i in range(1, N):
        newVert = mesh.addVertex(mesh.VPos[v1.ID, :]+direc*frac*i)
        verts.append(newVert)
    verts.append(v2)
    return verts

#Helper function for getBoxMesh
def addFaceTiles(mesh, stepSize, ebott, eright, etop, eleft):
    topRow = etop
    index = 1
    for index in range(1, len(eleft)):
        bottomRow = None
        if index == len(eleft)-1:
            bottomRow = ebott
        else:
            bottomRow = makeBoxEdge(mesh, eleft[index], eright[index], stepSize)
        #Now add the square faces on this part
        for i in range(0, len(topRow)-1):
            mesh.addFace([bottomRow[i], bottomRow[i+1], topRow[i+1], topRow[i]])
        topRow = bottomRow

#L is length along z
#W is width along x
#H is height along y
#stepSize is the length of each square tile.  By default there are no tiles
#(stepSize = -1).  If one of the sides is not an integer multiple of the step size,
#then round to the nearest step size that would make it an integer multiple along
#that dimension
def getBoxMesh(L = 1.0, W = 1.0, H = 1.0, C = np.array([0, 0, 0]), stepSize = -1):
    mesh = PolyMesh()
    endpoints = []
    for dZ in [L/2.0, -L/2.0]:
        for dH in [-H/2.0, H/2.0]:
            for dW in [-W/2.0, W/2.0]:
                endpoints.append(mesh.addVertex(C+np.array([dW, dH, dZ])))
    edgeIndices = [[0, 1], [1, 3], [3, 2], [2, 0], [1, 5], [5, 7], [7, 3], [7, 6], [6, 2], [0, 4], [4, 6], [4, 5]]
    edges = []
    edgesRev = []
    for edgePointers in edgeIndices:
        [v1, v2] = [endpoints[edgePointers[0]], endpoints[edgePointers[1]]]
        edges.append(makeBoxEdge(mesh, v1, v2, stepSize))
    for edge in edges:
        edgeRev = edge[:]
        edgeRev.reverse()
        edgesRev.append(edgeRev)
    #def addFaceTiles(mesh, stepSize, ebott, eright, etop, eleft):
    #Front Face
    addFaceTiles(mesh, stepSize, edges[0], edgesRev[1], edgesRev[2], edges[3])
    #Back Face
    addFaceTiles(mesh, stepSize, edgesRev[11], edgesRev[10], edges[7], edgesRev[5])
    #Left Face
    addFaceTiles(mesh, stepSize, edgesRev[9], edges[3], edges[8], edgesRev[10])
    #Right Face
    addFaceTiles(mesh, stepSize, edges[4], edgesRev[5], edgesRev[6], edgesRev[1])
    #Top Face
    addFaceTiles(mesh, stepSize, edgesRev[2], edges[6], edgesRev[7], edges[8])
    #Bottom Face
    addFaceTiles(mesh, stepSize, edges[11], edges[4], edges[0], edges[9])
    return mesh

def getRectMesh(P0, P1, P2, P3, stepSize = -1):
    mesh = PolyMesh()
    endpoints = [P0, P1, P2, P3]
    for i in range(0, len(endpoints)):
        endpoints[i] = mesh.addVertex(endpoints[i])
    edgeIndices = [[0, 1], [2, 1], [3, 2], [3, 0]]
    edges = []
    for edgePointers in edgeIndices:
        [v1, v2] = [endpoints[edgePointers[0]], endpoints[edgePointers[1]]]
        edges.append(makeBoxEdge(mesh, v1, v2, stepSize))
    addFaceTiles(mesh, stepSize, edges[0], edges[1], edges[2], edges[3])
    return mesh

def getTetrahedronMesh():
    mesh = PolyMesh()
    v1 = mesh.addVertex(np.array([-1, 1, 1]))
    v2 = mesh.addVertex(np.array([1, -1, 1]))
    v3 = mesh.addVertex(np.array([1, 1, -1]))
    v4 = mesh.addVertex(np.array([-1, -1, -1]))
    mesh.addFace([v1, v2, v3])
    mesh.addFace([v2, v4, v3])
    mesh.addFace([v3, v4, v1])
    mesh.addFace([v4, v2, v1])
    return mesh

def getOctahedronMesh():
    mesh = PolyMesh()
    v1 = mesh.addVertex(np.array([0, 0, 1]))
    v2 = mesh.addVertex(np.array([1, 0, 0]))
    v3 = mesh.addVertex(np.array([0, 1, 0]))
    v4 = mesh.addVertex(np.array([0, -1, 0]))
    v5 = mesh.addVertex(np.array([0, 0, -1]))
    v6 = mesh.addVertex(np.array([-1, 0, 0]))
    #Top Part
    mesh.addFace([v3, v1, v2])
    mesh.addFace([v3, v2, v5])
    mesh.addFace([v3, v5, v6])
    mesh.addFace([v3, v6, v1])
    #Bottom Part
    mesh.addFace([v1, v4, v2])
    mesh.addFace([v2, v4, v5])
    mesh.addFace([v5, v4, v6])
    mesh.addFace([v6, v4, v1])
    return mesh

def getIcosahedronMesh():
    mesh = PolyMesh()
    phi = (1+math.sqrt(5))/2
    #Use the unit cube to help construct the icosahedron
    #Front cube face vertices
    FL = mesh.addVertex(np.array([-0.5, 0, phi/2]))
    FR = mesh.addVertex(np.array([0.5, 0, phi/2]))
    #Back cube face vertices
    BL = mesh.addVertex(np.array([-0.5, 0, -phi/2]))
    BR = mesh.addVertex(np.array([0.5, 0, -phi/2]))
    #Top cube face vertices
    TF = mesh.addVertex(np.array([0, phi/2, 0.5]))
    TB = mesh.addVertex(np.array([0, phi/2, -0.5]))
    #Bottom cube face vertices
    BF = mesh.addVertex(np.array([0, -phi/2, 0.5]))
    BB = mesh.addVertex(np.array([0, -phi/2, -0.5]))
    #Left cube face vertices
    LT = mesh.addVertex(np.array([-phi/2, 0.5, 0]))
    LB = mesh.addVertex(np.array([-phi/2, -0.5, 0]))
    #Right cube face vertices
    RT = mesh.addVertex(np.array([phi/2, 0.5, 0]))
    RB = mesh.addVertex(np.array([phi/2, -0.5, 0]))
    
    #Add the icosahedron faces associated with each cube face
    #Front cube face faces
    mesh.addFace([TF, FL, FR])
    mesh.addFace([BF, FR, FL])
    #Back cube face faces
    mesh.addFace([TB, BR, BL])
    mesh.addFace([BB, BL, BR])
    #Top cube face faces
    mesh.addFace([TB, TF, RT])
    mesh.addFace([TF, TB, LT])
    #Bottom cube face faces
    mesh.addFace([BF, BB, RB])
    mesh.addFace([BB, BF, LB])
    #Left cube face faces
    mesh.addFace([LB, LT, BL])
    mesh.addFace([LT, LB, FL])
    #Right cube face faces
    mesh.addFace([RT, RB, BR])
    mesh.addFace([RB, RT, FR])
    
    #Add the icosahedron faces associated with each cube vertex
    #Front of cube
    mesh.addFace([FL, TF, LT]) #Top left corner
    mesh.addFace([BF, LB, FL]) #Bottom left corner
    mesh.addFace([FR, RT, TF]) #Top right corner
    mesh.addFace([BF, RB, FR]) #Bottom right corner
    #Back of cube
    mesh.addFace([LT, TB, BL]) #Top left corner
    mesh.addFace([BL, LB, BB]) #Bottom left corner
    mesh.addFace([RT, BR, TB]) #Top right corner
    mesh.addFace([BB, RB, BR]) #Bottom right corner
    
    return mesh

def getDodecahedronMesh():
    #Use icosahedron dual to help construct this
    icosa = getIcosahedronMesh()
    mesh = PolyMesh()
    #Add the vertex associated with each icosahedron face
    for f in icosa.faces:
        f.V = mesh.addVertex(f.getCentroid())
    #Add the face associated with each icosahedron vertex
    for v in icosa.vertices:
        verts = [f.V for f in v.getAttachedFaces()]
        vertsP = [mesh.VPos[V.ID, :] for V in verts]
        comparator = PointsCCWComparator(np.array([0, 0, 0]), vertsP[0])
        vertsP = [(i, vertsP[i]) for i in range(len(vertsP))]
        #Sort vertices in CCW order
        verts = [verts[x[0]] for x in sorted(vertsP, cmp=comparator.compare)]
        mesh.addFace(verts)
    return mesh

def getHemiOctahedronMesh():
    mesh = PolyMesh()
    v1 = mesh.addVertex(np.array([0, 0, 1]))
    v2 = mesh.addVertex(np.array([1, 0, 0]))
    v3 = mesh.addVertex(np.array([0, 1, 0]))
    v4 = mesh.addVertex(np.array([0, -1, 0]))
    v6 = mesh.addVertex(np.array([-1, 0, 0]))
    #Top Part
    mesh.addFace([v3, v1, v2])
    mesh.addFace([v3, v6, v1])
    #Bottom Part
    mesh.addFace([v1, v4, v2])
    mesh.addFace([v6, v4, v1])
    return mesh

def getSphereMesh(R, nIters):
    mesh = getOctahedronMesh()
    for i in range(nIters):
        mesh.evenTriangleRemesh()
        #Move points so that they're R away from the origin
        mesh.VPos = R*mesh.VPos/np.reshape(np.sqrt(np.sum(mesh.VPos**2, 1)), [mesh.VPos.shape[0], 1])
    return mesh

def getHemiSphereMesh(R, nIters):
    mesh = getHemiOctahedronMesh()
    for i in range(nIters):
        mesh.evenTriangleRemesh()
        #Move points so that they're R away from the origin
        mesh.VPos = mesh.VPos/np.reshape(np.sqrt(np.sum(mesh.VPos**2, 1)), [mesh.VPos.shape[0], 1])
    return mesh    

if __name__ == '__main__':
    icosahedronMesh = getIcosahedronMesh()
    icosahedronMesh.saveOffFile('icosahedron.off')
    dodecahedronMesh = getDodecahedronMesh()
    dodecahedronMesh.saveOffFile('dodecahedron.off')
    sphereMesh = getSphereMesh(1, 2)
    sphereMesh.saveOffFile("sphere.off")
    boxMesh = getBoxMesh(1, 1, 1, np.array([0, 0, 0]), 1)
    boxMesh.saveOffFile("box.off")
