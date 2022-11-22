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

        self.vx = np.ones((self.NX, self.NY)) * 5
        self.vx_n = np.ones((self.NX, self.NY))

        self.vy = np.ones((self.NX, self.NY)) * 2.5
        self.vy_n = np.ones((self.NX, self.NY))

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
        self.diffuse(self.dye, self.dye_n)
        self.dye[:] = self.dye_n[:]
        self.dissolve_dye()
        self.set_dye_bnd()
    
    def diffuse(self, y, y_n):
        a = self.diff * self.dt * self.cellSizeX * self.cellSizeY
        for k in range(20):
            for i in range(1, self.NX-1):
                for j in range(1, self.NY-1):
                    y_n[i, j] = (y[i, j] \
                        + a * 0.25 * (y_n[i+1, j]
                        + y_n[i-1, j]
                        + y_n[i, j+1]
                        + y_n[i, j-1])) / (1 + a)

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

        self.div[1:self.NX-1, 1:self.NY-1] = 0
    
    def set_within_bnd(self, idx, n):
        if idx > n-1:
            idx = n-1
        if idx < 0:
            idx = 0
        return idx

    def set_dye_bnd(self):
        self.dye[1:-1, 0] = self.dye[1:-1, 1]
        self.dye[1:-1, -1] = self.dye[1:-1, -2]
        self.dye[0, :] = self.dye[1, :]
        self.dye[-1, :] = self.dye[-2, :]
    
    def dissolve_dye(self):
        self.dye *= (1-self.dissolveRate)
        self.dye[self.dye<0.01] = 0

if __name__ == "__main__":
    cfd = Solver()
    cfd.run_solver()    