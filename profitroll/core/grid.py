import numpy as np

class Grid():
	""" This class encodes the 2D spatial grid used for the whole simulation.

	:param Lx: Length (in m) of the x-axis of the rectangle grid 
	:type Lx: float
	:param Ly: Length (in m) of the y-axis of the rectangle grid 
	:type Ly: float
	:param Nx: Number of subdivisions of the x-axis of the rectangle grid 
	:type Nx: int
	:param Ny: Number of subdivisions of the y-axis of the rectangle grid 
	:type Ny: int
	"""
	def __init__(self, Lx, Ly, Nx, Ny, **kwargs):
		""" Constructor method
		"""
		self.Lx = Lx
		self.Ly = Ly
		self.Nx = Nx
		self.Ny = Ny
        
		self.dx = self.Lx / self.Nx
		self.dy = self.Ly / self.Ny     
        
		grid = np.mgrid[0:self.Lx:self.dx, 0:self.Ly:self.dy]
		self.x_grid = grid[0,:,:]
		self.y_grid = grid[1,:,:]
