#TODO: Fix EPS weirdness
EPS = 1e-12
EPS_AREPLANAR = 1e-5
M_PI = 3.1415925
import math
import numpy as np
import numpy.linalg as linalg
import matplotlib.pyplot as plt


#############################################################
####                 PRIMITIVE CLASSES                  #####
#############################################################

class Plane3D(object):
    #P0 is some point on the plane, N is the normal
    def __init__(self, P0, N):
        self.P0 = np.array(P0)
        self.N = normalizeVec(N)
        self.resetEquation()

    def resetEquation(self):
        self.D = -self.P0.dot(self.N)

    def initFromEquation(self, A, B, C, D):
        N = np.array([A, B, C])
        self.P0 = (-D/N.dot(N))*np.array([A, B, C])
        self.N = normalizeVec(N)
        self.resetEquation()

    def distFromPlane(self, P):
        return self.N.dot(P) + self.D

    def __str__(self):
        return "Plane3D: %g*x + %g*y + %g*z + %g = 0"%(self.N[0], self.N[1], self.N[2], self.D)

class Line3D(object):
    def __init__(self, P0, V):
        self.P0 = np.array(P0)
        self.V = np.array(V)

    def intersectPlane(self, plane):
        P0 = plane.P0
        N = plane.N
        P = self.P0
        V = self.V
        if abs(N.dot(V)) < EPS:
            return None
        t = (P0.dot(N) - N.dot(P)) / (N.dot(V))
        intersectP = P + t*V
        return [t, intersectP]
    
    def intersectOtherLineRet_t(self, other):
        #Solve for (s, t) in the equation P0 + t*V0 = P1+s*V1
        #This is three equations (x, y, z components) in 2 variables (s, t)
        #User cramer's rule and the fact that there is a linear
        #dependence that only leaves two independent equations
        #(add the last two equations together)
        #[a b][t] = [e]
        #[c d][s]    [f]
        P0 = self.P0
        V0 = self.V
        P1 = other.P0
        V1 = other.V
        a = V0[0]+V0[2]
        b = -(V1[0]+V1[2])
        c = V0[1] + V0[2]
        d = -(V1[1]+V1[2])
        e = P1[0] + P1[2] - (P0[0] + P0[2])
        f = P1[1] + P1[2] - (P0[1] + P0[2])
        #print "[%g %g][t] = [%g]\n[%g %g][s]   [%g]"%(a, b, e, c, d, f)
        detDenom = a*d - c*b
        #Lines are parallel or skew
        if abs(detDenom) < EPS:
            return None
        detNumt = e*d - b*f
        detNums = a*f - c*e
        t = float(detNumt) / float(detDenom)
        s = float(detNums) / float(detDenom)
        #print "s = %g, t = %g"%(s, t)
        return (t, P0 + t*V0)
    
    def intersectOtherLine(self, other):
        ret = self.intersectOtherLineRet_t(other)
        if ret:
            return ret[1]
        return None
    
    def __str__(self):
        return "Line3D: %s + t%s"%(self.P0, self.V)


class Ray3D(object):
    def __init__(self, P0, V):
        self.P0 = np.array(P0)
        self.V = normalizeVec(V)
        self.line = Line3D(self.P0, self.V)
    
    def Copy(self):
        return Ray3D(self.P0, self.V)
    
    def Transform(self, matrix):
        self.P0 = mulHomogenous(matrix, self.P0.flatten())
        self.V = matrix[0:3, 0:3].dot(self.V)
        self.V = normalizeVec(self.V)
    
    def intersectPlane(self, plane):
        intersection = self.line.intersectPlane(plane)
        if intersection:
            if intersection[0] < 0:
                return None
            return intersection
    
    def intersectMeshFace(self, face):
        facePlane = face.getPlane()
        intersection = self.intersectPlane(facePlane)
        if not intersection:
            return None
        [t, intersectP] = intersection
        #Now check to see if the intersection is within the polygon
        #Do this by verifying that intersectP is on the same side
        #of each segment of the polygon
        verts = face.getVerticesPos()
        if verts.shape[0] < 3:
            return None
        lastCross = np.cross(verts[1, :]-verts[0, :], intersectP - verts[1, :])
        lastCross = normalizeVec(lastCross)
        for i in range(1, verts.shape[0]):
            v0 = verts[i, :]
            v1 = verts[(i+1)%verts.shape[0]]
            cross = np.cross(v1 - v0, intersectP - v1)
            cross = normalizeVec(cross)
            if cross.dot(lastCross) < EPS: #The intersection point is on the outside of the polygon
                return None
            lastCross = cross
        return [t, intersectP]

    def __str__(self):
        return "Ray3D: %s + t%s"%(self.P0, self.V)
    
#Axis-aligned bounding box class
class BBox3D(object):
    def __init__(self):
        self.b = np.array([[np.inf, np.inf, np.inf], [-np.inf, -np.inf, -np.inf]])
    
    def getDiagLength(self):
        dB = self.b[1, :] - self.b[0, :]
        return np.sqrt(dB.dot(dB))
    
    def getCenter(self):
        return np.mean(self.b, 0)
    
    def addPoint(self, P):
        self.b[0, :] = np.min((P, self.b[0, :]), 0)
        self.b[1, :] = np.max((P, self.b[1, :]), 0)
    
    def fromPoints(self, Ps):
        self.b[0, :] = np.min(Ps, 0)
        self.b[1, :] = np.max(Ps, 0)
    
    def Union(self, other):
        self.b[0, :] = np.min(self.b[0, :], other.b[0, :])
        self.b[1, :] = np.max(self.b[1, :], other.b[1, :])
    
    def __str__(self):
        coords = self.b.T.flatten()
        ranges = (self.b[1, :] - self.b[0, :]).flatten()
        return "BBox3D: [%g, %g] x [%g, %g] x [%g, %g],  Range (%g x %g x %g)"%tuple(coords.tolist() + ranges.tolist())

#############################################################
####                UTILITY FUNCTIONS                   #####
#############################################################
def splitIntoRGBA(val):
    A = (0xff000000&val)>>24
    R = (0x00ff0000&val)>>16
    G = (0x0000ff00&val)>>8
    B = (0x000000ff&val)
    return [R, G, B, A]

def extractFromRGBA(R, G, B, A):
    return ((A<<24)&0xff000000) | ((R<<16)&0x00ff0000) | ((G<<8)&0x0000ff00) | (B&0x000000ff)

def normalizeVec(V):
    return V/np.sqrt(np.sum(V**2))

#P is Nx3 matrix with N Points, one per row
def mulHomogenous(M, P):
    if len(P.shape) == 1:
        P = np.reshape(P, [1, 3])
    N = P.shape[0]
    PH = np.concatenate((P, np.ones((N, 1))), 1)
    ret = (M.dot(PH.T)).T
    return ret[:, 0:3]

#Project numpy vector V onto W
def projVec(V, W):
    return (V.dot(W)/(W.dot(W)))*W

def angleBetween(V, W):
    cosA = V.dot(W)/np.sqrt(V.dot(V)*W.dot(W))
    return np.arccos(cosA)

#Return the cosine of the angle between P1 and P2 with respect
#to "Vertex" as their common, shared vertex
def COSBetween(Vertex, P1, P2):
	V1 = P1 - Vertex
	V2 = P2 - Vertex
	dot = V1.dot(V2)
	magProduct = math.sqrt(np.sum(V1*V2)*np.sum(V2*V2))
	if (magProduct < EPS):
		return 0
	return float(dot) / float(magProduct)

class PointsCCWComparator(object):
    def __init__(self, C, VFirst):
        self.C = C #Center of reference for comparison
        self.VFirst = VFirst #First vertex in polygon
        self.N = VFirst - C #Normal direction
    
    def compare(self, V1, V2):
        a = V2[1] - self.VFirst
        b = V1[1] - self.VFirst
        triNormal = np.cross(a, b)
        dot = np.dot(triNormal, self.N)
        if dot > 0:
            return 1
        elif dot == 0:
            return 0
        return -1

#CCW rotate a vector (or point) V by angle "theta" around a line through P0 
#whose direction is specified by axis
def rotateAroundAxis(P0, axis, theta, V):
    #print "P0 = %s, axis = %s, V = %s"%(P0, axis, V)
    diffV = V - P0
    parV = projVec(diffV, axis) #Part of v along axis unaffected by rotation
    perpV = diffV - parV
    if perpV.dot(perpV) < EPS: #Hardly any perpendicular component
        return V
    u = normalizeVec(perpV)
    v = normalizeVec(axis)
    w = np.cross(u, v)
    
    #Put perpV into a frame where the rotation is about (0, 1, 0)
    fromFrame = np.zeros((3, 3))
    fromFrame[:, 0] = u
    fromFrame[:, 1] = v
    fromFrame[:, 2] = w
    toFrame = fromFrame.T
    perpV = toFrame.dot(perpV)
    #Rotate perpV by theta in that frame
    (cosTheta, sinTheta) = (np.cos(theta), np.sin(theta))
    rotMatrix = np.array([[cosTheta, 0, sinTheta], [0, 1, 0], [-sinTheta, 0, cosTheta]])
    perpV = rotMatrix.dot(perpV)
    #Bring perpV back into world coordinates and compute result
    perpV = fromFrame.dot(perpV)
    return P0 + perpV + parV

#Return True if the vertices in "verts" all lie
#in the same plane and False otherwise
def arePlanar(verts):
    if len(verts) <= 3:
        return True
    v0 = verts[1, :] - verts[0, :]
    v1 = verts[2, :] - verts[0, :]
    n = normalizeVec(np.cross(v0, v1))
    for i in range(3, verts.shape[0]):
        v = verts[i, :] - verts[0, :]
        v = normalizeVec(v)
        if n.dot(n) == 0:
            #If the first few points happened to be colinear
            n = v0.cross(v)
            if n.dot(n) == 0:
                #Still colinear...
                continue
        if abs(v.dot(n)) > EPS_AREPLANAR:
            return False
    return True

#If the vertices in "verts" form a convex 2D polygon 
#(in the order specified) return True.  Return False otherwise
def are2DConvex(verts):
    if len(verts) <= 3:
        return True
    if not arePlanar(verts):
        return False
    v0 = verts[0, :]
    v1 = verts[1, :]
    v2 = verts[2, :]
    lastCross = np.cross(v1-v0, v2-v1)
    for i in range(3, verts.shape[0]+1):
        v0 = v1
        v1 = v2
        v2 = verts[i%verts.shape[0]]
        cross = np.cross(v1-v0, v2-v1)
        if cross.dot(lastCross) < 0:
            return False
        lastCross = cross
    return True

#General purpose method for returning the normal of a face
#Assumes "verts" are planar and not all collinear
def getFaceNormal(verts):
    #This properly handles the case where three vertices
    #are collinear right after one another
    for i in range(2, verts.shape[0]):
        v1 = verts[i-1, :] - verts[0, :]
        v2 = verts[i, :] - verts[0, :]
        ret = np.cross(v1, v2)
        v1L2 = v1.dot(v1)
        v2L2 = v2.dot(v2)
        retL2 = ret.dot(ret)
        if v1L2 > 0 and v2L2 > 0 and retL2/(v1L2*v2L2) > EPS:
            return normalizeVec(ret)
    return None

#This function assumes the polygon is convex
def getPolygonArea(verts):
    if len(verts) < 3:
        return 0.0
    v1 = verts[1, :] - verts[0, :]
    v2 = verts[1, :] - verts[0, :]
    area = 0.0
    #Triangulate and add area of each triangle
    for i in range(2, len(verts)):
        v1 = v2
        v2 = verts[i, :] - verts[0, :]
        area = area + 0.5*np.sqrt(np.sum(np.cross(v1, v2)**2))
    return area

#Points of triangle are in rows of triangle
def getTriCircumcenter(tri):
    dV1 = tri[0, :] - tri[2, :]
    dV2 = tri[1, :] - tri[2, :]
    N = np.cross(dV1, dV2)
    P0 = tri[2, :] + 0.5*dV1
    P1 = tri[2, :] + 0.5*dV2
    V0 = np.cross(dV1, N)
    V1 = np.cross(dV2, N)
    a = V0[0]+V0[2]
    b = -(V1[0]+V1[2])
    c = V0[1] + V0[2]
    d = -(V1[1]+V1[2])
    e = P1[0] + P1[2] - (P0[0] + P0[2])
    f = P1[1] + P1[2] - (P0[1] + P0[2])
    #print "[%g %g][t] = [%g]\n[%g %g][s]   [%g]"%(a, b, e, c, d, f)
    detDenom = a*d - c*b
    #Lines are parallel or skew
    if abs(detDenom) < EPS:
        return None
    detNumt = e*d - b*f
    t = float(detNumt) / float(detDenom)
    return P0 + t*V0

#Get the closest point on the face represented by Vs to P
def getClosestPoint(Vs, P):
    #TODO: Finish translating this to numpy
    if len(Vs) < 3:
        return
    N = getFaceNormal(Vs)
    dV = P - Vs[0, :]
    #First find the closest point on the plane to this point
    PPlane = Vs[0, :] + N.projPerp(P)
    #Now check to see if PPlane is on the interior of the polygon
    ccws = np.zeros((len(Vs), 1))
    for i in range(len(Vs)):
        A = Vs[i]
        B = Vs[(i+1)%len(Vs)]
        v1 = PPlane - A
        v2 = PPlane - B
        M = np.ones((3, 3))
        M[1, :] = np.array([v1[0], v1[1], v1[2]])
        M[2, :] = np.array([v2[0], v2[1], v2[2]])
        ccw = linalg.det(M)
        if ccw < 0:
            ccws[i] = -1
        elif ccw > 0:
            ccws[i] = 1
        else:
            ccws[i] = 0
    if abs(sum(ccws)) == (ccws != 0).sum():
        return PPlane
    #Otherwise, the point is on the outside of the polygon
    #Check every edge to find the closest point
    minDist = np.inf
    ret = PPlane
    for i in range(len(Vs)):
        A = Vs[i]
        B = Vs[(i+1)%len(Vs)]
        dL = B - A
        PProjRel = dL.proj(PPlane - A)
        PClosest = A + PProjRel
        #Parameter representing where the point is
        #in the interior of the segment AB
        t = dL.Dot(PProjRel)/dL.squaredMag()
        #If the projected point is to the left of A, A is the closest
        if t < 0:
            PClosest = A
        #If the projected point is to the right of B, B is the closest
        elif t > 1:
            PClosest = B
        #Otherwise it is in the interior
        distSqr = (P - PClosest).squaredMag()

        if distSqr < minDist:
            minDist = distSqr
            ret = PClosest
    return ret

def PointsEqual(A, B):
    return abs(A[0] - B[0]) < EPS and abs(A[1] - B[1]) < EPS and abs(A[2] - B[2]) < EPS

if __name__ == '__main__':
    tri = np.zeros((3, 3))
    tri[:, 0:2] = np.random.randn(3, 2)
    print tri
    c = getTriCircumcenter(tri)
    print c
    plt.scatter(tri[:, 0], tri[:, 1], 20, 'b')
    plt.hold(True)
    for k in range(3):
        idx = [k, (k+1)%3]
        plt.plot(tri[idx, 0], tri[idx, 1], 'k')
    plt.scatter(c[0], c[1], 40, 'r')
    dR = tri - c
    R = np.sqrt(np.sum(dR**2, 1))
    print R
    tsine = np.linspace(0, 2*np.pi, 100)
    plt.plot(c[0] + R[0]*np.cos(tsine), c[1] + R[0]*np.sin(tsine), 'b')
    plt.show()

if __name__ == '__main__2':
    print "LINE INTERSECTION TEST"
    P0 = np.array([-2.5, 0, -2.5])
    V0 = np.array([-0, -0, -1])
    P1 = np.array([-2.5, -2.5, 0])
    V1 = np.array([0, -1, 0])
    line1 = Line3D(P0, V0)
    line2 = Line3D(P1, V1)
    intersection = line1.intersectOtherLine(line2)
    print intersection
    
    print "AXIS ROTATION TEST"
    np.random.seed(100)
    [R, S, V] = np.linalg.svd(np.random.randn(3, 3))
    PointRot = np.array([0, 0, 0])
    AxisRot = np.array([1, 0, 0])
    P = np.array([1, 2, 3])
    PRot1 = rotateAroundAxis(PointRot, AxisRot, 0.5, P)
    P = R.dot(P)
    AxisRot = R.dot(AxisRot)
    PRot2 = rotateAroundAxis(PointRot, AxisRot, 0.5, P)
    PRot2 = (R.T).dot(PRot2)
    print "PRot1 = ", PRot1
    print "PRot2 = ", PRot2
    
    P = np.array([1, 0, 0])
    P = rotateAroundAxis(np.array([0, 0, 0]), np.array([0, 0, 1]), np.pi/4, np.array([1, 0, 0]))
    print P
    
    
    print "AXIS ROTATION + LINE SEGMENT INTERSECTION TEST"
    P0 = np.array([1, 4, 0])
    P1 = np.array([5, 2, 0])
    P2 = np.array([2, 0, 0])
    P3 = np.array([3, 5, 0])
    PointRot = np.array([0, 0, 0])
    AxisRot = np.array([1, 1, 1])
    Angle = 0.5
    P0 = rotateAroundAxis(PointRot, AxisRot, Angle, P0)
    P1 = rotateAroundAxis(PointRot, AxisRot, Angle, P1)
    P2 = rotateAroundAxis(PointRot, AxisRot, Angle, P2)
    P3 = rotateAroundAxis(PointRot, AxisRot, Angle, P3)
    print "P0 = ",P0
    print "P1 = ",P1
    print "P2 = ",P2
    print "P3 = ",P3
    V0 = P1 - P0
    V1 = P3 - P2
    line1 = Line3D(P0, V0)
    line2 = Line3D(P2, V1)
    intersection = line1.intersectOtherLine(line2)
    intersection = rotateAroundAxis(PointRot, AxisRot, -Angle, intersection)
    print intersection
    
    print "CONVEXITY TEST"
    P = Plane3D(np.array([1, 1, 1]), np.array([1, 2, 3]))
    print P
    angle = 30
    angle = angle*3.141/180.0
    (cosA, sinA) = [math.cos(angle), math.sin(angle)]
    A = np.array([[cosA, -sinA, 0, 0],
                [sinA, cosA, 0, 0],
                [0, 0, 1, 0], [0, 0, 0, 1]])
    verts = np.array([[0, 0, 0], [0, 1, 0], [1, 1, 0], [2, 0.5, 0], [3, 0, 0], [1, -1, 0]])
    for i in range(verts.shape[0]):
        verts[i, :] = rotateAroundAxis(np.array([1, 1, 0]), np.array([1, 0, 0]), angle, verts[i, :])
    for v in verts:
        print v
    print are2DConvex(verts)
    
    print "BBOX Test"
    b = BBox3D(np.array([[1, 2, 3], [1, 2, 3]]))
    b.addPoint(np.array([0, 0, 0]))
    b.addPoint([-1, 5, -10])
    print b
