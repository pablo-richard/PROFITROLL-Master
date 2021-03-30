import numpy as np
from netCDF4 import Dataset

def create_initial_netcdf(path, Lx, Ly, Nx, Ny, dt, nb_state):
    """Creates an initial NetCDF for test cases of tropopause intrusion evolution.

    :param path: path where the NetCDF file will be created
    :type path: str
    :param Lx: Horizontal length of the grid
    :type Lx: float
    :param Ly: Vertical length of the grid
    :type Ly: float
    :param Nx: Horizontal number of cells of the grid
    :type Nx: int
    :param Ny: Vertical number of cells of the grid
    :type Ny: int
    :return: The created dataset
    :rtype: Dataset at NETCDF4 format
    """

#CREATION OF THE NETCDF FILE --------------------------------------------------
    handle = Dataset(path, 'w',format='NETCDF4')

#DIMENSIONS -------------------------------------------------------------------
    handle.createDimension("Nx", Nx)
    handle.createDimension("Ny", Ny)
    handle.createDimension("Nt", nb_state) 
    handle.createDimension("one", 1)

#ATTRIBUTE --------------------------------------------------------------------   
    # Geometry
    handle.Lx = Lx
    handle.Ly = Ly
    handle.Nx = Nx
    handle.Ny = Ny
    handle.dx = Lx / Nx
    handle.dy = Ly / Ny 
    
    ## Parameters
    handle.z_star = -500
    handle.gamma_1 = -4e-3 #K.m-1
    handle.gamma_2 = -8.5e-3 #K.m-1
    handle.Delta_zc = 500 
    handle.Delta_Tc = -5
    handle.g = 9.81
    handle.N_t = 0.01 # Brunt-Vaisala frequency of the troposphere (s^{-1})
    handle.N_s = 2e-2 # Brunt-Vaisala frequency of the stratosphere (s^{-1})
    handle.theta_00 = 300

#VARIABLES --------------------------------------------------------------------
    # "f8" is a data type: 64-bit floating point variable
    handle.createVariable("ut", "f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("vt", "f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("us", "f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("vs", "f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("w", "f8", ("Nx", "Ny", "Nt"))

    handle.createVariable("theta_t", "f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("Delta_T_bb", "f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("Delta_T_hist", "f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("Delta_z", "f8", ("Nx", "Ny", "Nt"))
    
    handle.createVariable("alpha_ut", "f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("alpha_vt", "f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("alpha_us", "f8", ("Nx", "Ny", "Nt"))
    handle.createVariable("alpha_vs", "f8", ("Nx", "Ny", "Nt"))
    
    handle.createVariable("t", "f8", ("Nt"))
    handle.createVariable("x_grid", "f8", ("Nx", "Ny"))
    handle.createVariable("y_grid", "f8", ("Nx", "Ny"))

#GEOMETRY INITIALIZATION ------------------------------------------------------
    grid = np.mgrid[0:handle.Lx:handle.dx, 0:handle.Ly:handle.dy]
    handle['x_grid'][:,:] = grid[0,:,:]
    handle['y_grid'][:,:] = grid[1,:,:]
    
#TIME INITIALIZATION-----------------------------------------------------------
    handle['t'][:] = dt * np.arange(nb_state)

    return handle


def v_stripe_test(path, Lx, Ly, Nx, Ny, dt, nb_state, dX, dY, depth=15):
    """ V stripe (intrusion) test case

    :param path: path where the NetCDF file will be created
    :type path: str
    :param Lx: Horizontal length of the grid
    :type Lx: float
    :param Ly: Vertical length of the grid
    :type Ly: float
    :param Nx: Horizontal number of cells of the grid
    :type Nx: int
    :param Ny: Vertical number of cells of the grid
    :type Ny: int
    :param dt: Time step
    :type dt: float
    :param nb_state: Number of created steps (same)
    :type nb_state: int
    :param dX: Nx - 2dX corresponds to the horizontal width of the v stripe
    :type dX: int
    :param dY: 2dY corresponds to the vertical width of the v stripe
    :type dY: int
    :param depth: Depth of the v stripe, default to 15
    :type depth: float, optional
    :return: The created dataset
    :rtype: Dataset at NETCDF4 format
    """
    
    handle = create_initial_netcdf(path, Lx, Ly, Nx, Ny, dt, nb_state)

#INITIALIZATION ---------------------------------------------------------------
        
    #Bubble creation
    [X,Y] = np.mgrid[0:Nx,0:Ny]
    F = np.zeros((Nx, Ny)) #handle.theta_00*np.ones((Nx, Ny))

    F[dX : Nx-dX, Ny//2 - dY : Ny//2 + dY] = \
        F[dX : Nx-dX, Ny//2 - dY : Ny//2 + dY] +\
        (np.abs(Y[dX : Nx-dX, Ny//2 - dY : Ny//2 + dY] - Ny//2)/dY -1)*depth
    
    left_ind= np.where(np.logical_and(\
                    (X - dX)**2 + (Y - Ny//2)**2 < (dY)**2,
                    X < dX ))
    F[left_ind] = F[left_ind] +\
             (np.sqrt((Y[left_ind]-Ny//2)**2 + (X[left_ind]-dX)**2 ) \
             /dY -1)*depth
                 
    right_ind= np.where(np.logical_and(\
                    (X - (Nx - dX - 1))**2 + (Y - Ny//2)**2 < dY**2,
                    X > (Nx - dX - 1) ))
    F[right_ind] = F[right_ind] +\
             (np.sqrt((Y[right_ind]-Ny//2)**2 + \
                      (X[right_ind]-(Nx-dX-1))**2 )/dY -1)*depth

    
    #Potential temperature
    handle['theta_t'][:,:,0] = F
    handle['theta_t'][:,:,1] = F
    
    handle['Delta_T_hist'][:,:,0] = handle.gamma_1 \
        * F * handle.g/(handle.N_t*handle.N_s*handle.theta_00)
    handle['Delta_T_hist'][:,:,1] =  handle['Delta_T_hist'][:,:,0] 
    handle['Delta_T_bb'][:,:,0] = handle['Delta_T_hist'][:,:,0]
    handle['Delta_z'][:,:,0] = np.zeros((Nx,Ny))
    handle['Delta_z'][:,:,1] = np.zeros((Nx,Ny))
    
    #Initial displacement guess for advection
    handle['alpha_ut'][:,:,0] = np.zeros((Nx,Ny))
    handle['alpha_vt'][:,:,0] = np.zeros((Nx,Ny))
    handle['alpha_us'][:,:,0] = np.zeros((Nx,Ny))
    handle['alpha_vs'][:,:,0] = np.zeros((Nx,Ny))
    
    #Uniform wind
    handle['ut'][:,:,0] = np.zeros((Nx,Ny))
    handle['vt'][:,:,0] = np.zeros((Nx,Ny))
    handle['us'][:,:,0] = np.zeros((Nx,Ny))
    handle['vs'][:,:,0] = np.zeros((Nx,Ny))
    handle['w'][:,:,0] = np.zeros((Nx,Ny))
    
    return handle

#------------------------------------------------------------------------------

def bubble_test(path, Lx, Ly, Nx, Ny, dt, nb_state, cx, cy, radius):
    """ Bubble test case (circular intrusion)

    :param path: path where the NetCDF file will be created
    :type path: str
    :param Lx: Horizontal length of the grid
    :type Lx: float
    :param Ly: Vertical length of the grid
    :type Ly: float
    :param Nx: Horizontal number of cells of the grid
    :type Nx: int
    :param Ny: Vertical number of cells of the grid
    :type Ny: int
    :param dt: Time step
    :type dt: float
    :param nb_state: Number of created steps (same)
    :type nb_state: int
    :param cx: Horizontal position of the center of the intrusion
    :type cx: int
    :param cy: Vertical position of the center of the intrusion
    :type cy: int
    :param radius: Radius of the circular intrusion
    :type radius: float
    :return: The created dataset
    :rtype: Dataset at NETCDF4 format
    """    
    handle = create_initial_netcdf(path, Lx, Ly, Nx, Ny, dt, nb_state)
        
    #Bubble creation
    bubble_indices = np.where ((handle['x_grid'][:,:] - cx)**2 +
                               (handle['y_grid'][:,:] - cy)**2 < radius**2 )
    
    F = np.zeros((Nx, Ny))
    F[bubble_indices] = 1
    
    #Potential temperature
    handle['theta_t'][:,:,0] = F
    handle['theta_t'][:,:,1] = F
    
    #Initial displacement guess for advection
    handle['alpha_u'][:,:,0] = np.zeros((Nx,Ny))
    handle['alpha_v'][:,:,0] = np.zeros((Nx,Ny))
    
    #Uniform wind
    handle['ut'][:,:,0] = np.zeros((Nx,Ny))
    handle['vt'][:,:,0] = np.zeros((Nx,Ny))

    return handle

#------------------------------------------------------------------------------

def gaussian_test(path, Lx, Ly, Nx, Ny, dt, nb_state):
    """ Gaussian (intrusion) test case of fixed variance Px=8 (horizontal) and Py=64 (vertical)

    :param path: path where the NetCDF file will be created
    :type path: str
    :param Lx: Horizontal length of the grid
    :type Lx: float
    :param Ly: Vertical length of the grid
    :type Ly: float
    :param Nx: Horizontal number of cells of the grid
    :type Nx: int
    :param Ny: Vertical number of cells of the grid
    :type Ny: int
    :param dt: Time step
    :type dt: float
    :param nb_state: Number of created steps (same)
    :type nb_state: int
    :return: The created dataset
    :rtype: Dataset at NETCDF4 format
    """  
    handle = create_initial_netcdf(path, Lx, Ly, Nx, Ny, dt, nb_state)
     
    #Gaussian creation
    thetatp = np.zeros((Nx, Ny))
    import scipy.ndimage as spnd 
    Px = 8 
    Py = 64
    thetaanom = 15
    for i in range(Nx): 
        for j in range(Ny):
            if abs(i-Nx/2) < Py and abs(j-Ny/2) < Px:
                thetatp[i,j] = thetaanom
    spnd.gaussian_filter(thetatp, 10, output=thetatp)
    
    #Potential temperature
    handle['theta_t'][:,:,0] = thetatp
    handle['theta_t'][:,:,1] = thetatp
    
    #Initial displacement guess for advection
    handle['alpha_ut'][:,:,0] = np.zeros((Nx, Ny))
    handle['alpha_vt'][:,:,0] = np.zeros((Nx, Ny))
    
    #Uniform wind
    handle['ut'][:,:,0] = np.zeros((Nx, Ny))
    handle['vt'][:,:,0] = np.zeros((Nx, Ny))
    

    return handle
    
#------------------------------------------------------------------------------
