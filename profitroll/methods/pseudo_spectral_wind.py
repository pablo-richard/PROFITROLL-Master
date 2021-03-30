from .spectral import geostwind

def pseudo_spectral_wind(history, grid, params, verbose, **kwargs):
    """Wrap the spectral methods to fit the architecture.
    
    :param history: Current history of state
    :type history: :class:`History` object
    :param grid: Spatial grid of the simulation
    :type grid: :class:`Grid` object:
    :param params: Dictionary of usefull parameters
    :type params: dictionary 
    :param verbose: verbose, defaults to 0
    :type verbose: int, optional
    """
    assert history.size > 0
    current_state = history.state_list[-1]
    
    ut, vt = geostwind(grid.Lx, grid.Ly, current_state.vrs['theta_t'], params, z=0, verbose=verbose)
    us, vs = geostwind(grid.Lx, grid.Ly, current_state.vrs['theta_t'], params, z=params['z_star'], verbose=verbose)
    
    current_state.vrs['ut'] = ut
    current_state.vrs['vt'] = vt
    current_state.vrs['us'] = us
    current_state.vrs['vs'] = vs