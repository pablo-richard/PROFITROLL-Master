import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import numpy as np
from netCDF4 import Dataset
from .animate import var2str

import ipywidgets as widgets

def plot_data(pathCDF, variable, time, cmap='magma'):
    """ Plot a variable from a NetCDF file using matlplotlib.

    :param pathCDF: path of the NetDCF file
    :type pathCDF: str
    :param variable: name of the variable to plot
    :type variable: str
    :param time: rank of the time step to plot
    :type time: int
    :param cmap: colormap to use, defaults to 'magma'
    :type cmap: str, optional
    """
    resultsCDF = Dataset(pathCDF, 'r', format='NETCDF4', parallel=False)
    
    # Get min and max value for the colorbar
    min_value = np.min(resultsCDF[variable][:])
    max_value = np.max(resultsCDF[variable][:])

    if(len(np.shape(resultsCDF[variable])) == 3):
    	# Figure Options
        fig = plt.figure(figsize=(12,8))
        ax = plt.subplot(111)

        ax.set_title(var2str(variable), fontsize=25, pad=15)

        # Color bar creation
        cbar = fig.colorbar(ScalarMappable(cmap=cmap, norm=Normalize(vmin=min_value, vmax=max_value)), ax=ax)
        cbar_ticks = cbar.get_ticks()
        cbar.set_ticks(cbar_ticks)

        # Axis tweaks
        ax.tick_params(labelsize=15, direction='in', length=10, width=1, pad=5, color='white')
        ax.set_ylabel(r'y Axis $(km)$', fontsize=15)
        ax.set_xlabel(r'x Axis $(km)$', fontsize=15)

    	# Subtitle
        htime = resultsCDF['t'][time]
        subtitle = 'Elapsed Time = {}'.format(int(np.floor(htime//3600))) + 'h {}min'.format(int(np.floor(htime%3600//60)))
        ax.text(0.1, -0.15, subtitle, size=plt.rcParams["axes.titlesize"], ha="center", transform=ax.transAxes)
        # Figure
        plt.imshow(resultsCDF[variable][:,:,time].T,
                    origin='lower', 
                    cmap=cmap,
                    vmin=min_value,
                    vmax=max_value)

    else :
        print(variable + ' is not a variable evolving on the 2D-grid through time.')
            
    resultsCDF.close()
    
def interactive_plot(pathCDF):
    """ Interactive plot of all the variables of a NetCDF file.

    :param pathCDF: Path to the NetCDF file
    :type pathCDF: str
    :return: the interactive object for interactive plots
    :rtype: :class:`ipywidgets.widgets.interactive` object
    """
    resultsCDF = Dataset(pathCDF, 'r', format='NETCDF4', parallel=False)

    times = resultsCDF['t'][:].data
    times_ind = np.arange(len(times))

    disp_times = [ '{}'.format(int(np.floor(times[k]//3600))) + 'h {}min'.format(int(np.floor(times[k]%3600//60))) for k in times_ind]
    
    vars_opts = list(resultsCDF.variables)
    cmap_opts = ['magma','Greys','hot','viridis','plasma','inferno','cividis']
    time_opts = [(disp_times[k], k) for k in times_ind]
    
    resultsCDF.close()

    vars_widg = widgets.ToggleButtons(options=vars_opts, 
                                    value=vars_opts[0], 
                                    description='variable :',
                                    disabled=False)
    
    time_widg = widgets.SelectionSlider(options=time_opts, 
    									value=0, 
    									description='t = ',
    									continuous_update=False, 
    									disabled=False)

    cmap_widg = widgets.Select(options=cmap_opts,
                                value='magma',
                                description='colormap :',
                                disabled=False)
    
    inter = widgets.interactive(plot_data, 
                                pathCDF=widgets.fixed(pathCDF),
                                cmap=cmap_widg,
                                variable=vars_widg,
                                time=time_widg)
    return inter
