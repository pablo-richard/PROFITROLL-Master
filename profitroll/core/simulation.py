import numpy as np
from netCDF4 import Dataset
import os
import time
from datetime import datetime
from copy import deepcopy

from .history import History
from .grid import Grid
from .netcdf_creator import create_results_netcdf, results_netcdf_frombackup

forced_attributes = ['T','Nt','methods','methods_kwargs','save_rate','backup_rate']

class Simulation():
    """ This class encodes the simulation

    :param history: List of states useful at each iteration (evolving through simulation)
    :type history: :class:`History` objects 
    :param params: Dictionary of specific parameters useful to the simulation
    :type params: dictionary
    :param grid: The 2D spatial grid used for the simulation
    :type grid: :class:`Grid` object
    :param methods: list of the methods used at each iteration of the simulations
    :type methods: list of functions
    :param methods_kwargs: list of dictionaries containing the arguments useful to each method used at each iteration of the simulations
    :type methods_kwargs: list of dictionaries
    :param output_folder: Path to the folder where the result and backup netCDF files will be saved
    :type output_folder: str
    :param T: list of the time duration of the simulations that have been launched,  defaults to empty list
    :type T: list of float, optional
    :param Nt: list of the number of time step of the simulations that have been launched,  defaults to empty list
    :type Nt: list of int, optional 
    :param save_rate:  list of the save rates of the simulations that have been launched,  defaults to empty list
    :type save_rate: list of int, optional
    :param backup_rate: list of the backup rates of the simulations that have been launched, defaults to empty list
    :type backup_rate: list of int, optional
    :type output_folder: str
    :param name: name of the simulation (which will contained in the result and backup files names), defaults to creation time
    :type name: str, optional
    :param saved_variables: list of the names of the variables that will be saved (if None all the variables will be saved), defaults to None
    :type saved_variables: list of str, optional
    :param verbose: Amount of informations that will be printed when running the simulation, defaults to 0
    :type verbose: int, optional
    """
    
    def __init__(self, initialCDF, methods, methods_kwargs, output_folder, save_rate=[], backup_rate=[], T=[], Nt=[], verbose=0, saved_variables=None, name=None, frombackup=False, pre_resultCDF=None):
        """ Constructor method
        
        :param initialCDF: netCDF file from which the parameters of the simulation, the initial history and the grid will be copied
        :type initialCDF: Dataset at NETCDF4 format
        :param frombackup: True if the simulation is created from a backup file (then a previous result file must be given), defaults to False
        :type frombackup: bool, optional
        :param pre_resultCDF: Previous result file if the simulation is created from a backup file (frombackup must be True), defaults to None
        :type pre_resultCDF: Dataset at NETCDF4 format, optional
        """
        if ((not frombackup) and (pre_resultCDF is not None)):
            initialCDF.close()
            pre_resultCDF.close()
            raise Exception('frombackup is False and 2 netCDF files were given ')

        # store the initial data
        self.history = History.fromCDF(initialCDF)
        self.params = {at: initialCDF.__dict__[at] for at in initialCDF.__dict__ if at not in forced_attributes}
        self.grid = Grid(**self.params)
        
        self.T = T
        self.Nt = Nt
        
        # checkout the methods
        self.methods = methods
        self.methods_kwargs = methods_kwargs if methods_kwargs is not None else len(methods)*[None]

        # handling the output netCDF files (save & backup)
        try:
            os.mkdir(output_folder)
        except FileExistsError:
            pass
        self.output_folder = output_folder
        self.save_rate = save_rate
        self.backup_rate = backup_rate
        self.saved_variables = saved_variables
        
        date = datetime.now()
        self.name = name if name is not None else date.strftime("%Y_%m_%d_%H:%M:%S")

        result_path = output_folder + '/results_'+self.name+'.nc'
        backup_path = output_folder + '/backup_'+self.name+'.nc'

        create_results_netcdf(result_path, initialCDF, **self.__dict__)
        if (frombackup and (pre_resultCDF is not None)):
            results_netcdf_frombackup(result_path, initialCDF, pre_resultCDF, **self.__dict__)
        create_results_netcdf(backup_path, initialCDF, **self.__dict__)


        initialCDF.close()
        if (pre_resultCDF is not None):
            pre_resultCDF.close()

        # other parameters
        self.verbose = verbose

    @classmethod
    def frombackup(cls, backupCDF, methods, methods_kwargs, output_folder, resultCDF=None, name=None, saved_variables=None, verbose=1):
        """ Other constructor method which construct a :class:`Simulation` object from a backup netCDF file. Informations to end the last simulation
        launched will be printed.

        :param backupCDF: 
        :type backupCDF: Dataset at NETCDF4 format
        :param methods: list of the methods used at each iteration of the simulations
	    :type methods: list of functions
        :param methods_kwargs: list of dictionaries containing the arguments useful to each method used at each iteration of the simulations
        :type methods_kwargs: list of dictionaries
        :param output_folder: Path to the folder where the result and backup netCDF files will be saved
        :type output_folder: str
        :param resultCDF: Previous result file, defaults to None
        :type resultCDF: Dataset at NETCDF4 format, optional
        :param name: name of the simulation (which will contained in the result and backup files names), defaults to creation time
	    :type name: str, optional
        :param saved_variables: list of the names of the variables that will be saved (if None all the variables will be saved), defaults to None
	    :type saved_variables: list of str, optional
        :param verbose: Amount of informations that will be printed when running the simulation, defaults to 1
	    :type verbose: int, optional
        """
        date = datetime.now()
        name = name if name is not None else 'frombackup_' + date.strftime("%Y_%m_%d_%H:%M:%S")
        
        T = deepcopy(backupCDF.T)
        Nt = deepcopy(backupCDF.Nt)

        save_rate = deepcopy(backupCDF.save_rate)
        backup_rate = deepcopy(backupCDF.backup_rate)

        missing_T = sum(T)-(backupCDF['t'][0]-resultCDF['t'][0]) if resultCDF is not None else sum(T)-backupCDF['t'][0]
        able = 2
        if (resultCDF is None) :
            able = 1 if len(backupCDF['t'][:].data) > 1 else 0
            dt = backupCDF['t'][1].data - backupCDF['t'][0].data if len(backupCDF['t'][:].data) > 1 else 1
            nb_v_step = backupCDF['t'][0].data / dt
            missing_Nt = sum(Nt)-nb_v_step
        else :
            valid_save = (np.where(resultCDF['t'][:].data < backupCDF['t'][0].data))[0] 
            nb_v_save = len(valid_save) 
            save_sim_comp = np.sum(np.array(Nt[:-1])//save_rate[:-1])
            missing_Nt = int(Nt[-1] - (nb_v_save - save_sim_comp)*save_rate[-1])

        T[-1] -= missing_T
        Nt[-1] -= missing_Nt

        if(able):
            print('To complete previous simulation, you might run this new simulation with :'
                    + '\n T = {}'.format(missing_T) 
                    + '\n Nt = {}'.format(missing_Nt) 
                    + '\n save_rate = {}'.format(save_rate[-1])
                    + '\n backup_rate = {}'.format(backup_rate[-1])
                    + '\n first_run = True')
            if(able==1):
                print('As no resultCDF has been given, this values are not certain. They are valid if you are using constant step methods and if the simulation began at time t = 0.')
        else:
            print('To complete previous simulation, you might run this new simulation with :' + '\n T = {}'.format(missing_T))
            print('As no resultCDF has been given and the backup file length correspond to only 1 step, the corresponding Nt cannot be retrive')

        return cls(backupCDF, methods, methods_kwargs, output_folder, save_rate, backup_rate, T=T, Nt=Nt,
                    verbose=verbose, saved_variables=saved_variables, 
                    name=name, frombackup=True, pre_resultCDF=resultCDF)


    def run(self, T, Nt, save_rate, backup_rate, first_run=True):
        """ Run a simulation for given time.

        :param T: Duration of the simulation
        :type T: float
        :param Nt: Number of time steps
        :type Nt: int
        :param save_rate: Rate of save in the result file
        :type save_rate: int
        :param backup_rate: Rate of refresh of the backup file
        :type backup_rate: int
        :param first_run: True if it is the first run of this simulation (for save purpose), defaults to True 
        :type first_run: bool, optional
        """
        if (backup_rate%save_rate):
            raise Exception('For recovery from backup purpose save_rate must divide backup_rate')

        cpu_tot_time = np.zeros(len(self.methods))
        simu_time = time.time()

        # Saving parameters of the new run
        backupCDF = Dataset(self.output_folder + '/backup_'+self.name+'.nc', 'r+', format='NETCDF4', parallel=False)
        resultsCDF = Dataset(self.output_folder + '/results_'+self.name+'.nc', 'r+', format='NETCDF4', parallel=False)
        for ob in [self, backupCDF, resultsCDF]:
            ob.T = np.append(ob.T, T)
            ob.Nt = np.append(ob.Nt, Nt)
            ob.save_rate = np.append(ob.save_rate, save_rate)
            ob.backup_rate = np.append(ob.backup_rate, backup_rate)
        backupCDF.close()
        resultsCDF.close()  

        if self.verbose:
            print("          ------------------------")
            print("          |  RUNNING SIMULATION  |")
            print("          ------------------------")
        for iter_nb in range(Nt):
            print("\n\nIteration ", iter_nb, "...") if self.verbose else None
            # first handle saving
            if (iter_nb % self.backup_rate[-1] == 0) and not (iter_nb==0 and not first_run):
                backupCDF = Dataset(self.output_folder + '/backup_'+self.name+'.nc', 'r+', format='NETCDF4', parallel=False)
                self.history.save(backupCDF, backup=True)
                backupCDF.close()
                print("---> backup refreshed at iteration "+str(iter_nb)) if self.verbose else None
            if iter_nb % self.save_rate[-1] == 0 and not (iter_nb==0 and not first_run):
                resultsCDF = Dataset(self.output_folder + '/results_'+self.name+'.nc', 'r+', format='NETCDF4', parallel=False)
                self.history.save(resultsCDF, backup=False, saved_variables=self.saved_variables)  
                resultsCDF.close()  
                print("---> saved results of iteration "+str(iter_nb)) if self.verbose else None

            # then perform forward
            cpu_time = self.forward()
            cpu_tot_time += cpu_time    
        
        # Last save/backup
        backupCDF = Dataset(self.output_folder + '/backup_'+self.name+'.nc', 'r+', format='NETCDF4', parallel=False)
        self.history.save(backupCDF, backup=True)
        backupCDF.close()
        resultsCDF = Dataset(self.output_folder + '/results_'+self.name+'.nc', 'r+', format='NETCDF4', parallel=False)
        self.history.save(resultsCDF, backup=False, saved_variables=self.saved_variables)  
        resultsCDF.close() 

        # FINAL PRINT : Print Total and Mean CPU time per method
        for ind, method in enumerate(self.methods):
            print("\n\nTotal CPU time for method ", method.__name__, " = {:.2f}".format(cpu_tot_time[ind]), " seconds") if self.verbose else None
            print("Mean CPU time for method ", method.__name__, " per call = {:.2f}".format(cpu_tot_time[ind]/self.Nt[-1]), " seconds") if self.verbose else None

        simu_time = time.time() - simu_time
        print("\n**************************************************\n")
        print("TOTAL METHODS TIME = {:.2f}".format(np.sum(cpu_tot_time)), " seconds")
        print("TOTAL SIMULATION TIME = {:.2f}".format(simu_time), " seconds")

    def forward(self):
        """ One step of simulation : apply each method to the history.
        """
        cpu_time = np.zeros(len(self.methods))
        for ind, (method, kwargs) in enumerate(zip(self.methods, self.methods_kwargs)):
            t0 = time.time()
            print("      *** Proceeding to method: "+method.__name__) if self.verbose > 1 else None
            method(**self.__dict__, **kwargs)
            cpu_time[ind] = time.time() - t0 
            print("      *** CPU time = {:.2f}".format(cpu_time[ind]), " seconds") if self.verbose > 1 else None
        return cpu_time
