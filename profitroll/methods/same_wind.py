def same_wind(history, **kwargs):
    """Set the tropopause wind at their previous values for the new state
    
    :param history: Current history of state
    :type history: :class:`History` object
    """
    assert history.size > 1
    previous_state = history.state_list[-2]
    current_state = history.state_list[-1]

    update_var = ['ut','vt','us','vs']
    for var in update_var:
        current_state.vrs[var] = previous_state.vrs[var]
		
