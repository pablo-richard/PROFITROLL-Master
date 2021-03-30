def end_pop(history, **kwargs):
    """ Pop the first state of the history
    
    :param history: Current history of state
    :type param: :class:`History` object
    """
    assert history.size > 0
    history.pop(0)

		
