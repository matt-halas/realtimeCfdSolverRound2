from tkinter.tix import CELL
from types import CellType
import pygame
import sys

import numpy as np

from utilities import *

CELL_SIZE = 20
NX = 40
dx = 1
NY = 40
dy = 1

X_center = np.arange(0, NX, 1, dtype=int)
Y_center = np.arange(0, NY, 1, dtype=int)

cell_idx = np.zeros((NX*NY, 2))
for j in range(NY):
    for i in range(NX):
        cell_idx[two_to_one(i, j)] = [X_center[i], Y_center[j]]


VISC = 1
DIFF = 0.1
dt = 0.1

dye_c = np.zeros(NX*NY)
dye_n = np.zeros(NX*NY)

angle = np.pi/6 #Angle of uniform flow field in radians
magnitude = 5 #Magnitude of uniform flow field
vx = np.ones(NX*NY) * magnitude * np.cos(angle)
vy = np.ones(NX*NY) * magnitude * np.sin(angle)

vx0 = np.zeros(NX*NY)
vy0 = np.zeros(NX*NY)

pygame.init()
screen = pygame.display.set_mode((NX * CELL_SIZE, NY * CELL_SIZE))
pygame.display.set_caption("This time it'll work")

def check_events():
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            sys.exit()
        elif event.type == pygame.KEYDOWN:
            if event.key == pygame.K_q:
                sys.exit()

    if pygame.mouse.get_pressed()[0]:
        add_dye()

def add_dye(dye_speed=100):
    mouse_pos = pygame.mouse.get_pos()
    #mouse_pos[0] is x position - differs from physical position outlined in cell centers
    x_idx = mouse_pos[0] // CELL_SIZE
    y_idx = mouse_pos[1] // CELL_SIZE
    dye_c[two_to_one(x_idx, y_idx)] += dye_speed

def draw_dye():
    # Uses the dye_c value at the end of each time step to draw on the display
    for i in range(len(dye_c)):
        rect = pygame.Rect((i % NX) * CELL_SIZE,
            (i // NY) * CELL_SIZE,
            CELL_SIZE, CELL_SIZE)
        if dye_c[i] > 255:
            rect_color = (255, 255, 255)
        else:
            rect_color = (dye_c[i],dye_c[i],dye_c[i])
        pygame.draw.rect(screen, rect_color, rect)

def diffuse():
    dye_n[:] = dye_c[:]
    calc_diffuse()
    set_diff_bnd()
    dye_c[:] = dye_n[:]

def calc_diffuse(iterations=10):
    # Try in a few different ways - Std diffusion, Gauss Siedel relaxation with d_n=0 and d_n=d_c
    # Dye_n is not solving correctly, probably missed the denominator
    for k in range(iterations):
        for i in range(1, NX-1):
            for j in range(1, NY-1):
                dye_n[two_to_one(i, j)] = (dye_c[two_to_one(i, j)] \
                    + DIFF * 0.25 * (dye_n[two_to_one(i+1, j)]
                    + dye_n[two_to_one(i-1, j)]
                    + dye_n[two_to_one(i, j+1)]
                    + dye_n[two_to_one(i, j-1)])) \
                        / (1 + DIFF)

def set_diff_bnd():
    # In this state, the boundaries keep all the dye in the simulation
    for i in range(1, NX-1):
        dye_n[two_to_one(i, 0)] = dye_n[two_to_one(i, 1)]
        dye_n[two_to_one(i, NY-1)] = dye_n[two_to_one(i, NY-2)]
    for j in range(1, NY-1):
        dye_n[two_to_one(0, j)] = dye_n[two_to_one(1, j)]
        dye_n[two_to_one(NX-1, j)] = dye_n[two_to_one(NX-2, j)]
    dye_n[two_to_one(0, 0)] = np.mean([dye_n[two_to_one(1, 0)], dye_n[two_to_one(0, 1)]])
    dye_n[two_to_one(NX-1, 0)] = np.mean([dye_n[two_to_one(NX-2, 0)], dye_n[two_to_one(NX-1, 1)]])
    dye_n[two_to_one(0, NY-1)] = np.mean([dye_n[two_to_one(0, NY-2)], dye_n[two_to_one(1, NY-1)]])
    dye_n[two_to_one(NX-1, NY-1)] = np.mean([dye_n[two_to_one(NX-2, NY-1)], dye_n[two_to_one(NX-1, NY-2)]])

def advect():
    for i in range(1, NX-1):
        for j in range(1, NY-1):
            x_cell, y_cell = cell_idx[two_to_one(i, j)] #cell_idx is only the index of a given cell
            x_cell*=dx #Multiply by dx and dy to get the actual physical location to make it compatible
            y_cell*=dy #with dx and dy other than 1
            x_vel = vx[two_to_one(i,j)]
            y_vel = vy[two_to_one(i,j)]
            x_adv = x_cell - x_vel*dt
            y_adv = y_cell - y_vel*dt
            x_idx = np.floor(x_adv / dx, dtype=int)
            y_idx = np.floor(y_adv / dy, dtype=int)
            #dye_n[two_to_one(i, j)] = 

            
    

def runSolver():
    while True:
        check_events()
        diffuse()
        advect()
        draw_dye()
        pygame.display.flip()

if __name__ == "__main__":
    mouseDown = True
    runSolver()

#TODO: Add more comments
#TODO: Add boundary conditions - Diffusion done
#TODO: Add units with physical meaning
#TODO: Add velocity field