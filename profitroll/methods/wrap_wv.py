from .advection_step_3P import advection_step_3P
from .spectral import vertwind
from ..core.state import State #, variables
import numpy as np

def wrap_wv(history, grid, params, alpha_method, order_alpha, F_method, verbose=0, **kwargs):
    """Wrap the water vapor method to fit the architecture.
    
    :param history: Current history of state
    :type history: :class:`History` object
    :param grid: Spatial grid of the simulation
    :type grid: :class:`Grid` object
    :param params: Dictionary of usefull parameters
    :type params: dictionary 
    :param alpha_method: see :class:`advection_step_3P`
    :type alpha_method: str
    :param order_alpha: see :class:`advection_step_3P`
    :type order_alpha: int
    :param F_method: see :class:`advection_step_3P`
    :type F_method: str
    :param verbose: verbose, defaults to 0
    :type verbose: int, optional
    """
    assert history.size > 2
    pre_state = history.state_list[-3]
    cur_state = history.state_list[-2]
    new_state = history.state_list[-1]

    dt = cur_state.t - pre_state.t # constant step
    
    invar = np.array([pre_state.vrs['Delta_z'],pre_state.vrs['Delta_T_hist']]) 
    a_us, a_vs, outvar = advection_step_3P(pre_state.vrs['alpha_us'],
                                           pre_state.vrs['alpha_vs'],
                                           invar,
                                           dt,
                                           cur_state.vrs['us'],
                                           cur_state.vrs['vs'],
                                           grid.dx,
                                           grid.dy,
                                           alpha_method,
                                           order_alpha,
                                           F_method,
                                           verbose)
    print("      us vs done") if verbose > 2 else None
    
    new_dz = outvar[0]
    new_dT_hist = outvar[1] 
    
    #UPDATE OF W ---------------------------------------------------
    k_hour = int(3600/dt)
    if ((np.floor(cur_state.t/dt)-1)%k_hour==0 ):
        cur_w = vertwind(grid.Lx, grid.Ly, cur_state.vrs['theta_t'], pre_state.vrs['theta_t'], dt, params, z=params['z_star'])
        new_w = vertwind(grid.Lx, grid.Ly, new_state.vrs['theta_t'], cur_state.vrs['theta_t'], dt, params, z=params['z_star'])
        mean_w = (cur_w + new_w)/2.
        cur_state.vrs['Delta_z'] += k_hour * dt * mean_w
        new_dz += k_hour * dt * mean_w
        
    dT_disp = params['gamma_2'] * cur_state.vrs['Delta_z']
    dT_cloud = params['Delta_Tc'] * ( cur_state.vrs['Delta_z'] > params['Delta_zc'] ) 
    dT_bb = cur_state.vrs['Delta_T_hist'] + dT_disp + dT_cloud

    cur_state.vrs['alpha_us'] = a_us
    cur_state.vrs['alpha_vs'] = a_vs
    
    new_state.vrs['Delta_z'] = new_dz 
    new_state.vrs['Delta_T_hist'] = new_dT_hist
    
    cur_state.vrs['Delta_T_bb'] = dT_bb
