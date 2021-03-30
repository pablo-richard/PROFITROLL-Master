from .advection_step_3P import advection_step_3P
from ..core.state import State

def wrap_advection_step_3P(history, grid, params, alpha_method, order_alpha, F_method, verbose=0, **kwargs):
    """Wrap the :class:`advection_step_3P` method to fit the architecture

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
    assert history.size > 1
    pre_state = history.state_list[-2]
    cur_state = history.state_list[-1]

    dt = cur_state.t - pre_state.t # constant step

    new_state = State.copy(cur_state)
    new_state.t += dt              # constant step
    
    a_ut, a_vt, theta_new = advection_step_3P(pre_state.vrs['alpha_ut'],
                                              pre_state.vrs['alpha_vt'],
                                              pre_state.vrs['theta_t'],
                                              dt,
                                              cur_state.vrs['ut'],
                                              cur_state.vrs['vt'],
                                              grid.dx,
                                              grid.dy,
                                              alpha_method,
                                              order_alpha,
                                              F_method,
                                              verbose)
    print("      ut vt done") if verbose > 2 else None
    
    cur_state.vrs['alpha_ut'] = a_ut
    cur_state.vrs['alpha_vt'] = a_vt
    
    new_state.vrs['theta_t'] = theta_new
    
    history.append(new_state)
