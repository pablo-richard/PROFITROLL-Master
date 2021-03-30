import numpy as np

from .upstream_interp import upstream_interp

def advection_step_3P(alpha_u_minus, alpha_v_minus, field_minus,
                      dt, u, v, dx, dy,
                      alpha_method,
                      order_alpha,
                      F_method,
                      verbose=0):
    
    if alpha_method == 'damped_bicubic' or F_method == 'damped_bicubic':
         a = 0.5
         B = 4.0
         d0 = 3.25E-5
         d = 0.5 * np.sqrt(
             np.square( (np.roll(u,-1,0) - np.roll(u,1,0)) / (2 * dx)   - \
                        (np.roll(v,-1,1) - np.roll(v,1,1)) / (2 * dy) ) + \
             np.square( (np.roll(u,-1,1) - np.roll(u,1,1)) / (2 * dy)   + \
                        (np.roll(v,-1,0) - np.roll(v,1,0)) / (2 * dx) ) )
         d2 = (d.copy())/d0
         d2[np.where(d2<1)] = 1
         f = a * d * d2**B
         kappa = f * dt / (1 + f * dt)
         print("kappa: ", np.mean(kappa)," , ", np.min(kappa)," , ", np.max(kappa)) if verbose > 2 else None 
        
    # ITERATIVE ESTIMATION OF THE DISPLACEMENT-------------------------------
    # The displacement alpha is estimated iteratively by interpolating the wind 
    # in x -alpha_minus, where alpha is the previous estimate. At each time
    # step, a number of iterations order_alpha is used, and the iterative 
    # scheme is initialized with the estimate at the previous time step.
    for k in range(order_alpha):
        # Staniforth et Al. states that for the interpolation of the 
        # estimated displacement, linear interpolation is usually 
        # sufficient. Bicubic interpolation is more accurate but more 
        # expensive. Damped bicubic interpolation also combines the accuracy 
        # of bicubic interpolation with a relaxation based on the linear 
        # interpolation to limit small scale numerical noise.
        if k < order_alpha-1 and alpha_method !='linear':
            method = 'linear'
        else:
            method = alpha_method
            
        print("      advection_step_3P with alpha order "+str(k)+" and method "\
              +method) if verbose > 2 else None
        
        if method == 'damped_bicubic':
            [alpha_u, alpha_v] = (dt/dx)* ( 
                kappa * upstream_interp(alpha_u_minus, alpha_v_minus, 
                        np.array([u,v]), method='linear',verbose=verbose ) 
                + (1- kappa) * upstream_interp(alpha_u_minus, alpha_v_minus, 
                            np.array([u,v]), method='bicubic',verbose=verbose))
        else:
            [alpha_u, alpha_v] = (dt/dx)*upstream_interp(alpha_u_minus,
                                                     alpha_v_minus,
                                                     np.array([u,v]),
                                                     method=method,
                                                     verbose=verbose)
        alpha_u_minus = alpha_u
        alpha_v_minus = alpha_v
        
#------------------------------------------------------------------------
    
    # The field at time tk + dt is udpated by interpolating the field at 
    # time tk -dt at the locations x - 2* alpha.
    #field_plus = upstream_interp(2*alpha_u, 2*alpha_v, field_minus,
    #                             method=F_method, verbose=verbose)

    if F_method == 'damped_bicubic':
        Ia = upstream_interp(2*alpha_u, 2*alpha_v, field_minus,
                                        method='bicubic',verbose=verbose )
        Id = upstream_interp(2*alpha_u, 2*alpha_v, field_minus,
                                        method='linear',verbose=verbose )
        field_plus =  kappa * Id +(1- kappa)* Ia 

 
    else:
        field_plus = upstream_interp(2*alpha_u, 2*alpha_v, field_minus,
                                        method=F_method, verbose=verbose )
    
    return alpha_u, alpha_v, field_plus