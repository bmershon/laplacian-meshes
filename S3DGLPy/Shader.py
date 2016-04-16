from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GL.ARB.framebuffer_object import *
from OpenGL.GL.EXT.framebuffer_object import *
from OpenGL.GL.ARB.vertex_buffer_object import *
from OpenGL.GL.ARB.geometry_shader4 import *
from OpenGL.GL.EXT.geometry_shader4 import *

from Cameras3D import *

class Shader(object):
    def __init__(self, vertexFile, fragmentFile, geometryFile):
        fin = open(vertexFile, "r")
        self.vSource = fin.read()
        fin.close()
        fin = open(fragmentFile, "r")
        self.fSource = fin.read()
        fin.close()
        fin = open(geometryFile, "r")
        self.gSource = fin.read()
        fin.close()
        self.prog = None
        
        self.uPMatrix = None
        self.uMVMatrix = None
        self.uAmbientColor = None
        self.uLightPos = None
        self.uLightColor = None
    
    def compile(self):
        #Compile shaders
        (self.vs, self.fs, self.gs) = (0, 0, 0)
        self.vs = glCreateShader(GL_VERTEX_SHADER)
        self.fs = glCreateShader(GL_FRAGMENT_SHADER)
        self.gs = glCreateShader(GL_GEOMETRY_SHADER_EXT)
        glShaderSource(self.vs, self.vSource)
        glShaderSource(self.fs, self.fSource)
        glShaderSource(self.gs, self.gSource)
        
        glCompileShader(self.vs)
        log = glGetShaderInfoLog(self.vs)
        if log:
            print "Vertex Shader: ", log
        
        glCompileShader(self.fs)
        log = glGetShaderInfoLog(self.fs)
        if log:
            print "Fragment Shader: ", log
        
        glCompileShader(self.gs)
        log = glGetShaderInfoLog(self.gs)
        if log:
            print "Geometry Shader: ", log
        
        #Compile program
        self.prog = glCreateProgram()
        glAttachShader(self.prog, self.vs)
        glAttachShader(self.prog, self.fs)
        glAttachShader(self.prog, self.gs)
        glLinkProgram(self.prog)
        glUseProgram(self.prog)
        
        #Get uniform locations
        self.uPMatrix = getUniformLocation(self.prog, "uPMatrix");
        self.uMVMatrix = getUniformLocation(self.prog, "uMVMatrix");
        self.uAmbientColor = getUniformLocation(self.prog, "uAmbientColor");
        self.uLightPos = getUniformLocation(self.prog, "uLightPos")
        self.uLightColor = getUniformLocation(self.prog, "uLightColor")
    
    def use(self, camera, ambientColor = np.array([0.5, 0.5, 0.5]), lightColor = np.array([0.5, 0.5, 0.5])):
        p = camera.getPerspectiveMatrix()
        glUniformMatrix4fv(self.uPMatrix, 1, False, p)
        m = camera.getModelviewMatrix()
        glUniformMatrix4fv(self.uMVMatrix, 1, False, m)
        glUniform3fv(self.uAmbientColor, 1, ambientColor)
        #Make a "head lamp" by default
        glUniform3fv(self.uLightPos, 1, camera.eye)
        glUniform3fv(self.uLightColor, 1, lightColor)
