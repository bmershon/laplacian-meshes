from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.arrays import vbo
import numpy as np
import scipy.io as sio

def getColorPickingTexture():
    J = sio.loadmat('colors.mat')
    J = J['J']
    texId = glGenTextures(1)
    glBindTexture(GL_TEXTURE_2D, texId)
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
    glTexImage2D(GL_TEXTURE_2D, 0, 3, J.shape[0], J.shape[1], 0, GL_RGB, GL_UNSIGNED_BYTE, J)
    return texId

def drawColorPicker(width, height, texID):
    glDisable(GL_DEPTH_TEST)
    glDisable(GL_LIGHTING)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluOrtho2D(0, width, 0, height)
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
    
    #Draw spectrogram
    glEnable(GL_TEXTURE_2D)
    glBindTexture(GL_TEXTURE_2D, texID)
    glColor3f(1, 1, 1)
    glBegin(GL_QUADS)
    glTexCoord2f(1.0, 0.0); glVertex2i(200, 0)
    glTexCoord2f(1.0,1.0); glVertex2i(200, 200)
    glTexCoord2f(0.0, 1.0); glVertex2i(0, 200)
    glTexCoord2f(0.0,0.0); glVertex2i(0, 0)
    glEnd()
    glDisable(GL_TEXTURE_2D)
    
    glEnable(GL_LIGHTING)
    glEnable(GL_DEPTH_TEST)
