import numpy as np

def upstream_interp(alpha_x, alpha_y, F, method='linear', verbose=0, ho=0.15, **kwargs):
    """
    upstream_interp interpolates a multidimensionnal field F from a 2D grid to an 'upstream' unstructured mesh defined by the displacements alpha_x, alpha_y. \
    If F were a continuous field, we would have: F_int(x,y) = F(x-alpha, y-alpha)
        
    :param alpha_x: a two dimensional field of displacement along the first dimension
    :type alpha_x: ndarray
    :param alpha_y: a two dimensional field of displacement along the second dimension
    :type alpha_y: ndarray
    :param F: 2D-Field to be advected. \
        A vectorial field can be advected (the X and Y dimension should come second and third)
    :type F: ndarray
    :param method: method used for the interpolation: 'nearest', 'linear', or 'bicubic'. Defaults to 'linear'
    :type method: string, optional
    :param verbose: a scalar between 0 and 2. The higher the number, the more prints. Defaults to 0
    :type verbose: int, optional
    :param ho: Distance around grid points where additional diffusion is added with the 'damping' method. Defaults to 0.15.
    :type ho: int, optional
    :raises "Unknown method for interpolation": Invalid string as a method for interpolation
    :return: F_int: Advected field. 
    :rtype: ndarray
    """


    print("         upstream_interp called with method: ", method) if verbose > 2 else None

    if len(F.shape)==2:
        F = np.array([F])
        
    [dim,Nx,Ny] = F.shape
    F_int = np.zeros((dim,Nx,Ny))
    if method=='nearest':
         # Coordinates of the points of the upstream mesh
         [x_grid, y_grid] = np.mgrid[0:Nx, 0:Ny]
         x_grid = x_grid - alpha_x
         y_grid = y_grid - alpha_y
         
         # Apply periodic boundary coniditions. Note that nothing is done
         # for x,y<0 as the indices are already periodic in Python. 
         x_grid[np.where(x_grid > Nx-1)] -= (Nx-1)
         y_grid[np.where(y_grid > Ny-1)] -= (Ny-1)
         
         # Fetch the closest neighbor
         F_int[:,:,:] = F[np.round(x_grid), np.round(y_grid),:] ##### seems to generate bugs
         
    elif method=='linear':
        [X, Y] = np.mgrid[0:Nx,0:Ny] 
        
        Xi = np.mod(X - alpha_x, Nx - 1)
        Yi = np.mod(Y - alpha_y, Ny - 1)
        
        Xt = np.ceil(Xi).astype(int)
        Yt = np.ceil(Yi).astype(int)
        
        Xc = Xt - Xi
        Yc = Yt - Yi
        
        F_int[:, X, Y] = Xc * Yc * F[:, Xt - 1, Yt -1] \
                   + Xc * (1 - Yc) * F[:, Xt - 1, Yt] \
                   + (1 - Xc) * (1 - Yc) * F[:, Xt, Yt] \
                   + (1 - Xc) * Yc * F[:, Xt, Yt - 1] 
                   
                   
    elif method=='diffusive':
        [X, Y] = np.mgrid[0:Nx,0:Ny] 
        
        Xi = X - alpha_x
        Yi = Y - alpha_y
        
        Xt = np.ceil(Xi).astype(int)
        Yt = np.ceil(Yi).astype(int)
        
        Xc = Xt - Xi
        Yc = Yt - Yi
        
        Xt = np.mod(Xt, Nx)
        Yt = np.mod(Yt, Ny)
        
        Ft = F[:, Xt, Yt]
        
        F_int = Xc * Yc * np.roll(Ft, (1, 1), (1,2)) \
                   + Xc * (1 - Yc) * np.roll(Ft, 1, 1) \
                   + (1 - Xc) * (1 - Yc) * Ft \
                   + (1 - Xc) * Yc * np.roll(Ft, 1, 2) 
        

    elif method=='bicubic':
         # For a given point, bi-cubic interpolation fits the value of the 
        # four surrounding points as well as the slope at each of these 
        # points. The slope is evaluated using centered finite differences.
        # We can reduce this operation to a weight for each of the sixteen
        # surrounding points. We compute these weights for each grid cell.
        # Each of the 16 is thus stored in an array..
        
        #-----------------------------------------------------------------
        [X, Y] = np.mgrid[0:Nx,0:Ny] 
        
        Xi = X - alpha_x
        Yi = Y - alpha_y
        
        Xf = np.mod( np.floor(Xi).astype(int), Nx)
        Yf = np.mod( np.floor(Yi).astype(int), Ny)
        
        Xb = np.abs(Xi - np.floor(Xi))
        Yb = np.abs(Yi - np.floor(Yi))
        
        
        
        # We define 'rolled' arrays that appears multiple times 
        F00 = np.roll(F, (1,1), (1,2))
        F01 = np.roll(F, (1,0), (1,2))
        F02 = np.roll(F, (1,-1), (1,2))
        F03 = np.roll(F, (1,-2), (1,2))
        F10 = np.roll(F, (0,1), (1,2))
        F11 = np.roll(F, (0,0), (1,2))
        F12 = np.roll(F, (0,-1), (1,2))
        F13 = np.roll(F, (0,-2), (1,2))
        F20 = np.roll(F, (-1,1), (1,2))
        F21 = np.roll(F, (-1,0), (1,2))
        F22 = np.roll(F, (-1,-1), (1,2))
        F23 = np.roll(F, (-1,-2), (1,2))
        F30 = np.roll(F, (-2,1), (1,2))
        F31 = np.roll(F, (-2,0), (1,2))
        F32 = np.roll(F, (-2,-1), (1,2))
        F33 = np.roll(F, (-2,-2), (1,2))
        
        # Computation of the weights -------------------------------------
        A00 = F11
        
        A01 = -0.5 * F10 + 0.5 * F12
        
        A02 = F10 - 2.5 * F11 + 2 * F12 -  0.5 * F13
            
        A03 = -.5 * F10 + 1.5 * F11 - 1.5 * F12 + 0.5 * F13
            
        A10 = -.5 * F01 + 0.5 * F21
        
        A11 = 0.25 * F00 - 0.25 * F02 - 0.25 * F20 + 0.25 * F22
            
        A12 = -.5 * F00 + 1.25 * F01 - F02 +  0.25 * F03 + 0.5 * F20 - 1.25  \
            * F21 + F22 - 0.25 * F23
                    
        A13 = 0.25 * F00 - 0.75 * F01 + 0.75 * F02 - 0.25 * F03 - 0.25 * F20 \
            + 0.75 * F21 - 0.75 * F22 + 0.25 * F23
                
        A20 = F01 - 2.5 * F11 + 2 * F21 -  0.5 * F31
            
        A21 = -.5 * F00 + 0.5 * F02 + 1.25 * F10 - 1.25 * F12 - F20 + F22 + \
            0.25 * F30 - 0.25 * F32
                    
        A22 = F00 - 2.5 * F01 + 2 * F02 - 0.5 * F03 - 2.5 * F10 +  6.25 * \
            F11 - 5 * F12 + 1.25 * F13 + 2 * F20  - 5 * F21 + 4 * F22 - \
            F23 - 0.5 * F30 + 1.25 * F31  - F32 + 0.25 * F33
                            
        A23 = -.5 * F00 + 1.5 * F01 - 1.5 * F02 +  0.5 * F03 + 1.25 * F10 - \
            3.75 * F11 + 3.75 * F12 -  1.25 * F13 - F20 + 3 * F21 - 3 * F22 \
                +  F23 + 0.25 * F30 - 0.75 * F31 + 0.75 * F32 - 0.25 * F33
                    
        A30 = -.5 * F01 + 1.5 * F11 - 1.5 * F21 + 0.5 * F31
            
        A31 = 0.25 * F00 - 0.25 * F02 - 0.75 * F10 + 0.75 * F12 + 0.75 * F20 \
            - 0.75 * F22 - 0.25 * F30 + 0.25 * F32;
        
        A32 = -.5 * F00 + 1.25 * F01 - F02 + 0.25 * F03 + 1.5 * F10 - 3.75 *\
            F11 + 3 * F12 - 0.75 * F13 - 1.5 * F20 + 3.75 * F21 - 3 * F22 +\
            0.75 * F23 + 0.5 * F30 - 1.25 * F31 + F32 - 0.25 * F33
		
        A33 = 0.25 * F00 - 0.75 * F01 + 0.75 * F02 - 0.25 * F03 - 0.75 * F10 \
            + 2.25 * F11 - 2.25 * F12 + 0.75 * F13 + 0.75 * F20 - 2.25 * F21 +\
            2.25 * F22 - 0.75 * F23 - 0.25 * F30 + 0.75 * F31 - 0.75 * F32 +\
            0.25 * F33
    
        
        #Update of F
        F_int = (A00[:,Xf, Yf] + A01[:,Xf, Yf] * Yb + A02[:,Xf, Yf] * Yb**2 + \
                 A03[:,Xf, Yf] * Yb**3) + \
                (A10[:,Xf, Yf] + A11[:,Xf, Yf] * Yb + A12[:,Xf, Yf] * Yb**2 + \
                 A13[:,Xf, Yf] * Yb**3) * Xb + \
                (A20[:,Xf, Yf] + A21[:,Xf, Yf] * Yb + A22[:,Xf, Yf] * Yb**2 + \
                 A23[:,Xf, Yf] * Yb**3) * Xb**2 + \
                (A30[:,Xf, Yf] + A31[:,Xf, Yf] * Yb + A32[:,Xf, Yf] * Yb**2 + \
                 A33[:,Xf, Yf] * Yb**3) * Xb**3  
    else:
        raise Exception("Unknown method for interpolation: " + method)
    if dim==1:
        F_int = F_int[0,:,:]
    
    return F_int

