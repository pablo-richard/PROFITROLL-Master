from .state import State

class History():
	""" This class encodes the list of :class:`State` objects useful at each simulation step (evolving through simulation)

	:param state_list: List of :class:`State` objects 
	:type state_list: List of :class:`State` objects 
	:param size: Length of state_list
	:type size: int
	"""
	def __init__(self, state_list):
		""" Constructor method

		:param state_list: List of :class:`State` objects 
		"""
		self.state_list = state_list
		self.size = len(state_list)
		
	@classmethod
	def fromCDF(cls, netCDF_file):
		""" Other constructor method which construct a :class:`History` object from a netCDF file

		:param netCDF_file: NetCDF file used to create the :class:`History` object
		:type netCDF_file: Dataset at NETCDF4 format
		"""
		size = netCDF_file.dimensions['Nt'].size
		if size > 0:
			state_list = [State.fromCDF(netCDF_file, k) for k in range(size)]
			return cls(state_list)
		else:
			raise 'Empty CDF while initialising History'

	def save(self, netCDF_file, backup=False, saved_variables=None):
		""" Saves the first state of state_list in a netCDF file (or the entire state_list if backup is True)

		:param netCDF_file: File where the state will be saved
		:type netCDF_file: Dataset at NETCDF4 format
		:param backup: if False, the first state of state_list is added to the netCDF file . if True the netCDF file state informations are replaced by the complete state_list, defaults to False
		:type backup: bool, optional
		:param saved_variables: List of the variables (str) which will be saved, if None, all the variables are saved, defaults to None
		:type saved_variables: list of str ,optional
		"""
		if not backup:
        # the oldest state is supposed to be completely known
        # thus it is this one we choose to export
			self.state_list[0].save(netCDF_file, backup=False, saved_vrs=saved_variables)
		else:
			for ind, state in enumerate(self.state_list):
				state.save(netCDF_file, backup=True, k=ind)
				
	def append(self, state):
		""" Adds a :class:`State` object at the end of the state_list

		:param state: This state will be added at the end of the state_list
		:type state: :class:`State` object
		"""
		self.state_list.append(state)
		self.size += 1
		
	def pop(self, k=None):
		""" Removes a :class:`State` object from the state_list
		
		:param k: The (k+1)-th element of the state_list will be removed, if None the first element is removed, defaults to None
		:type k: int, optional
		"""
		try:
			self.state_list.pop(k)
			self.size -= 1
		except Exception as error:
			raise error
