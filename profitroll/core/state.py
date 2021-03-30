import numpy as np
from copy import deepcopy
        
forced_variables = ['x_grid','y_grid','t'] # General variables     

class State():
    """ This class encodes an object which contains all the variables at a given time t.

	:param t: Corresponding time
	:type t: float
	:param vrs: Dictionary of the variables and their values
	:type vrs: dictionary
	"""
    def __init__(self, t, vrs={}):
        """Constructor method
        """
        self.t = t
        self.vrs = deepcopy(vrs)

    @classmethod
    def fromCDF(cls, netCDF_file, k=None):
        """ Other constructor method which construct a :class:`State` object from a netCDF file

        :param netCDF_file: NetCDF file used to create the :class:`State` object
		:type netCDF_file: Dataset at NETCDF4 format
        :param k: time rank of the state to create in the NetCDF file, defaults to None
        :type k: int, optional
        """
        t = np.float(netCDF_file['t'][k].data) if k is not None else None

        variables = [var for var in netCDF_file.variables if var not in forced_variables]
        vrs = {var: netCDF_file[var][:,:,k].data if k is not None else netCDF_file[var][:].data
               for var in variables}

        return cls(t, vrs)

    @classmethod
    def copy(cls, otherState):
        """Creates a copy of an other :class:`State` object

        :param otherState: :class:`State` object to copy
		:type otherState: :class:`State` object
        """
        return cls(otherState.t, otherState.vrs)

    def save(self, netCDF_file, saved_vrs=None, backup=False, k=None):
        """Save the state into a given NetCDF file.

        :param netCDF_file: NetCDF file where the :class:`State` object will be saved
		:type netCDF_file: Dataset at NETCDF4 format
        :param saved_vrs: list of the variables to save, defaults to None (all the variables are saved)
		:type saved_vrs: list of str, optional
        :param backup: True if this is the NetCDF backup file (the previous states are erased), else the state is simply addded, defaults to False
        :type backup: bool, optional
        :param k: time rank given to the state in the netCDF file , defaults to None (added at the end)
        :type k: int, optional
        """
        # current location to fill
        k = netCDF_file.dimensions['Nt'].size if not backup else k

        # filling the variables
        netCDF_file['t'][k] = self.t

        variables = [var for var in netCDF_file.variables if var not in forced_variables]
        export_variables = saved_vrs if saved_vrs is not None else variables
        for var in export_variables:
            netCDF_file[var][:,:,k] = self.vrs[var]
