class Settings:
    # Class that stores the settings for the solver and physical
    # properties such as time step, diffusion, viscosity
    def __init__(self):
        self.screenHeight = 800
        self.screenWidth = 800
        self.NX = 80
        self.NY = 80
        self.cellSizeX = self.screenWidth / self.NX
        self.cellSizeY = self.screenHeight / self.NY