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

cell_center = np.zeros((NX*NY, 2))
for j in range(NY):
    for i in range(NX):
        cell_center[two_to_one(i, j)] = [X_center[i], Y_center[j]]
cell_idx = np.zeros((NX*NY, 2))
for j in range(NY):
    for i in range(NX):
        cell_idx[two_to_one(i, j)] = [X_center[i], Y_center[j]]

diffusivity = 1
viscosity = 1
dt = 0.1

dye_c = np.zeros(NX*NY)
dye_n = np.zeros(NX*NY)

vx_n = np.ones(NX*NY)
vy_n = np.ones(NX*NY)

vx_c = np.zeros(NX*NY)
vy_c = np.zeros(NX*NY)

p = np.zeros(NX*NY)
div = np.zeros(NX*NY)

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
    vx_c[two_to_one(x_idx, y_idx)] = 2
    vy_c[two_to_one(x_idx, y_idx)] = 2

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

def diffuse(y_c, y_n, diff):
    # Try in a few different ways - Std diffusion, Gauss Siedel relaxation with d_n=0 and d_n=d_c
    # Dye_n is not solving correctly, probably missed the denominator
    a = diff * dt * dx * dy
    for k in range(20):
        for i in range(1, NX-1):
            for j in range(1, NY-1):
                y_n[two_to_one(i, j)] = (y_c[two_to_one(i, j)] \
                    + a * 0.25 * (y_n[two_to_one(i+1, j)]
                    + y_n[two_to_one(i-1, j)]
                    + y_n[two_to_one(i, j+1)]
                    + y_n[two_to_one(i, j-1)])) \
                        / (1 + a)
    #Current diffusion routine was too slow, going back to the bad one
'''    for i in range(1, NX-1):
        for j in range(1, NY-1):
            y_n[two_to_one(i, j)] = y_c[two_to_one(i, j)] \
                + a * 0.25 * (y_c[two_to_one(i+1, j)]
                + y_c[two_to_one(i-1, j)]
                + y_c[two_to_one(i, j+1)]
                + y_c[two_to_one(i, j-1)]
                - 4 * y_c[two_to_one(i, j)])'''

def advect(y_c, y_n):
    for i in range(1, NX-1):
        for j in range(1, NY-1):
            x_cell, y_cell = cell_center[two_to_one(i, j)]
            x_vel = vx_c[two_to_one(i,j)]
            y_vel = vy_c[two_to_one(i,j)]
            x_adv = x_cell - x_vel*dt
            y_adv = y_cell - y_vel*dt
            xi = np.int(np.floor(x_adv / dx))
            yi = np.int(np.floor(y_adv / dy))
            interpx1 = lerp(cell_center[two_to_one(xi, yi), 0],
                        cell_center[two_to_one(xi+1, yi), 0],
                        x_adv, y_c[two_to_one(xi, yi)],
                        y_c[two_to_one(xi+1, yi)])
            interpx2 = lerp(cell_center[two_to_one(xi, yi+1), 0],
                        cell_center[two_to_one(xi+1, yi+1), 0],
                        x_adv, y_c[two_to_one(xi, yi+1)],
                        y_c[two_to_one(xi+1, yi+1)])
            y_n[two_to_one(i, j)] = lerp(cell_center[two_to_one(xi, yi), 1],
                                           cell_center[two_to_one(xi, yi+1), 1],
                                           y_adv, interpx1, interpx2)

def project(vx_c, vy_c):
    p[:] = 0
    for i in range(1, NX-1):
        for j in range(1, NY-1):
            div[two_to_one(i, j)] = -(vx_c[two_to_one(i+1, j)] - vx_c[two_to_one(i-1, j)]) / (2 * dx) \
                - (vy_c[two_to_one(i, j+1)] - vy_c[two_to_one(i, j-1)]) / (2 * dy)
    
    for k in range(20):
        for i in range(1, NX-1):
            for j in range(1, NY-1):
                p[two_to_one(i, j)] = (div[two_to_one(i, j)] + p[two_to_one(i-1, j)] + p[two_to_one(i+1, j)] \
                    + p[two_to_one(i, j+1)] + p[two_to_one(i, j-1)]) / 4
    
    for i in range(1, NX-1):
        for j in range(1, NY-1):
            vx_c -= (p[two_to_one(i+1, j)] - p[two_to_one(i-1, j)]) / (2 * dx)
            vy_c -= (p[two_to_one(i, j+1)] - p[two_to_one(i, j-1)]) / (2 * dy)

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
    dye_c[:] = dye_n[:]

def set_vel_bnd():
    for i in range(0, NX):
        vx_c[two_to_one(i, 0)] = 0
        vx_c[two_to_one(i, NY-1)] = 0
        vy_c[two_to_one(i, 0)] = 0
        vy_c[two_to_one(i, NY-1)] = 0
    for j in range(0, NY):
        vx_c[two_to_one(0, j)] = 0
        vx_c[two_to_one(NX-1, j)] = 0
        vy_c[two_to_one(0, j)] = 0
        vy_c[two_to_one(NX-1, j)] = 0

def stepDye():
    dye_n[:] = dye_c[:]
    diffuse(dye_c, dye_n, diffusivity)
    dye_c[:] = dye_n[:]
    advect(dye_c, dye_n)
    dye_c[:] = dye_n[:]
    set_diff_bnd()

def stepVel():
    set_vel_bnd()
    vx_n[:], vy_n[:] = [vx_c[:], vy_c[:]]
    diffuse(vx_c, vx_n, viscosity)
    diffuse(vy_c, vy_n, viscosity)
    vx_c[:], vy_c[:] = [vx_n[:], vy_n[:]]
    advect(vx_c, vx_n)
    advect(vy_c, vy_n)
    vx_c[:], vy_c[:] = [vx_n[:], vy_n[:]]
    project(vx_c, vy_c)
    set_vel_bnd()
    

def runSolver():
    while True:
        check_events()
        stepDye()
        stepVel()
        draw_dye()
        pygame.display.flip()

if __name__ == "__main__":
    mouseDown = True
    runSolver()

#TODO: Add more comments
#TODO: Add boundary conditions - Diffusion done
#TODO: Add units with physical meaning
