import pygame
import numpy as np
import sys

from utilities import *
from settings import Settings

class Solver:

    def __init__(self):
        pygame.init()
        self.settings = Settings()
        self.NX = self.settings.NX
        self.NY = self.settings.NY
        self.cellSizeX = self.settings.cellSizeX
        self.cellSizeY = self.settings.cellSizeY

        self.cellCenter = np.zeros((NX, NY, 2))
        self.populate_cell_center()

        self.diff = 0.001
        self.dt = 0.1
        self.dissolveRate = 0.001

        self.dye = np.zeros((self.NX, self.NY))
        self.dye_n = np.zeros((self.NX, self.NY))

        self.vx = np.zeros((self.NX, self.NY))
        self.vx_n = np.zeros((self.NX, self.NY))

        self.vy = np.zeros((self.NX, self.NY))
        self.vy_n = np.zeros((self.NX, self.NY))

        self.p = np.zeros((self.NX, self.NY))
        self.div = np.zeros((self.NX, self.NY))

        self.screen = pygame.display.set_mode((self.settings.screenWidth,
            self.settings.screenHeight))
        pygame.display.set_caption("This time it'll work right")
    
    def populate_cell_center(self):
        X_center = np.linspace(self.cellSizeX/2,
            self.NX*self.cellSizeX-self.cellSizeX/2,
            self.NX)
        Y_center = np.linspace(self.cellSizeY/2,
            self.NY*self.cellSizeY-self.cellSizeY/2,
            self.NY)
        for j in range(self.NY):
            for i in range(self.NX):
                self.cellCenter[i, j] = [X_center[i], Y_center[j]]

    def run_solver(self):
        while True:
            self.check_events()
            self.step_dye()
            self.step_vel()
            self.update_screen()

    def check_events(self):
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                sys.exit()
            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_q:
                    sys.exit()
        
        if pygame.mouse.get_pressed()[0]:
            self.add_dye()
    
    def add_dye(self, dyeSpeed = 200):
        mouse_pos = pygame.mouse.get_pos()
        #mouse_pos[0] is x position - differs from physical position outlined in cell centers
        x_idx = np.int(mouse_pos[0] // self.cellSizeX)
        y_idx = np.int(mouse_pos[1] // self.cellSizeY)
        self.dye[x_idx, y_idx] += dyeSpeed

    def update_screen(self):
        #Uses the dye value at the end of each time step to draw on the display
        for i in range(self.NX):
            for j in range(self.NY):
                rect = pygame.Rect(i * self.cellSizeX, j * self.cellSizeY,
                    self.cellSizeX, self.cellSizeY)
                if self.dye[i, j] > 255:
                    rect_color = (255, 255, 255)
                else:
                    rect_color = (self.dye[i, j],
                        self.dye[i, j],self.dye[i, j])
                pygame.draw.rect(self.screen, rect_color, rect)
        
        pygame.display.flip()

    def step_dye(self):
        self.dye_n[:] = self.dye[:]
        self.advect(self.dye, self.dye_n)
        self.dye[:] = self.dye_n[:]
        self.set_cont_bnd(self.dye)
        self.diffuse(self.dye, self.dye_n)
        self.dye[:] = self.dye_n[:]
        self.dissolve_dye()
        self.set_cont_bnd(self.dye)
    
    def step_vel(self):
        self.vx_n[:] = self.vx[:]
        self.vy_n[:] = self.vy[:]
        self.diffuse(self.vx, self.vx_n, isVel=True)
        self.diffuse(self.vy, self.vy_n, isVel=True)
        self.set_vel_bnd()
        self.project()
        self.vx_n[:] = self.vx[:]
        self.vy_n[:] = self.vy[:]
        self.advect(self.vx, self.vx_n)
        self.advect(self.vy, self.vy_n)
        self.project()
        self.set_vel_bnd()
    
    def diffuse(self, y, y_n, isVel=False):
        a = self.diff * self.dt * self.cellSizeX * self.cellSizeY
        for k in range(10):
            for i in range(1, self.NX-1):
                for j in range(1, self.NY-1):
                    y_n[i, j] = (y[i, j] \
                        + a * 0.25 * (y_n[i+1, j]
                        + y_n[i-1, j]
                        + y_n[i, j+1]
                        + y_n[i, j-1])) / (1 + a)
            if isVel:
                self.set_vel_bnd()
            else:
                self.set_cont_bnd(y_n)

    def advect(self, y, y_n):
        for i in range(1, self.NX-1):
            for j in range(1, self.NY-1):
                x_cell, y_cell = self.cellCenter[i,j]
                x_adv = x_cell - self.vx[i,j]*self.dt
                y_adv = y_cell - self.vy[i,j]*self.dt
                xi = np.int(np.floor(x_adv/self.cellSizeX))
                yi = np.int(np.floor(y_adv/self.cellSizeY))
                xi = self.set_within_bnd(xi, self.NX)
                yi = self.set_within_bnd(yi, self.NY)
                interpx1 = lerp(self.cellCenter[xi, yi, 0],
                    self.cellCenter[xi+1, yi, 0],
                    x_adv, y[xi, yi], y[xi+1, yi])
                interpx2 = lerp(self.cellCenter[xi, yi+1, 0],
                    self.cellCenter[xi+1, yi+1, 0],
                    x_adv, y[xi, yi+1], y[xi+1,yi+1])
                y_n[i,j] = lerp(self.cellCenter[xi, yi, 1],
                    self.cellCenter[xi, yi+1, 1], y_adv,
                    interpx1, interpx2)
    
    def project(self):
        self.p[:] = 0
        self.div[:] = 0

        self.div[1:self.NX-1, 1:self.NY-1] = (self.vx[2:, 1:self.NY-1] \
            - self.vx[0:-2, 1:self.NY-1]) / (2 * self.cellSizeX) \
            + (self.vy[1:self.NX-1, 2:] - self.vy[1:self.NX-1, 0:-2]) \
            / (2 * self.cellSizeY)
        
        self.set_cont_bnd(self.div)
        self.set_cont_bnd(self.p)
        
        for k in range(10):
            for i in range(1, NX-1):
                for j in range(1, NY-1):
                    self.p[i,j] = (self.div[i,j] + self.p[i+1,j] + self.p[i-1,j]
                        + self.p[i,j+1] + self.p[i,j-1]) / 4
            self.set_cont_bnd(self.p)
        
        self.vx[1:self.NX-1, 1:self.NY] -= (self.p[2:, 1:self.NY]
            - self.p[:-2, 1:self.NY]) / 2
        self.vy[1:self.NX-1, 1:self.NY] -= (self.p[1:self.NX, 2:]
            - self.p[1:self.NX, :-2]) / 2
        self.set_vel_bnd()
    
    def set_within_bnd(self, idx, n):
        if idx > n-1:
            idx = n-1
        if idx < 0:
            idx = 0
        return idx

    def set_cont_bnd(self, x):
        x[1:-1, 0] = x[1:-1, 1]
        x[1:-1, -1] = x[1:-1, -2]
        x[0, 1:-1] = x[1, 1:-1]
        x[-1, 1:-1] = x[-2, 1:-1]
        self.set_corner_bnd(x)
    
    def set_vel_bnd(self):
        self.vx[1:-1, 0] = self.vx[1:-1, 1]
        self.vx[1:-1, -1] = self.vx[1:-1, -2]
        self.vx[0, 1:-1] = -self.vx[1, 1:-1]
        self.vx[-1, 1:-1] = -self.vx[-2, 1:-1]
        self.vy[1:-1, 0] = -self.vy[1:-1, 1]
        self.vy[1:-1, -1] = -self.vy[1:-1, -2]
        self.vy[0, 1:-1] = self.vy[1, 1:-1]
        self.vy[-1, 1:-1] = self.vy[-2, 1:-1]
        self.set_corner_bnd(self.vx)
        self.set_corner_bnd(self.vy)
    
    def set_corner_bnd(self, x):
        x[0, 0] = (x[1, 0] + x[0, 1]) / 2
        x[-1, 0] = (x[-2, 0] + x[-1, 1]) / 2
        x[0, -1] = (x[0, -2] + x[1, -1]) / 2
        x[-1, -1] = (x[-2, -1] + x[-1, -2]) / 2
    
    def dissolve_dye(self):
        self.dye *= (1-self.dissolveRate)
        self.dye[self.dye<0.01] = 0

if __name__ == "__main__":
    cfd = Solver()
    cfd.run_solver()    