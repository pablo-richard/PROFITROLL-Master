import matplotlib.pyplot as plt
from matplotlib.animation import ArtistAnimation
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import numpy as np
from netCDF4 import Dataset
from IPython.display import HTML
from base64 import b64encode
import os

def make_video(pathCDF, save_path, variable, cmap='magma'):
    """Loads the data stored in a NetCDF file, builds a video of the asked variable and saves it. It returns an HTML animation that can be displayed in a Jupyter notebook.
    
    :param pathCDF: path of the NetCDF file where to read the data
    :type pathCDF: string
    :param save_path: path where to save the video
    :type save_path: string
    :param variable: physical quantity to plot on video, must belong to the NetCDF file
    :type variable: string
    :return: an HTML animation containing the video asked
    :rtype: IPython.core.display.HTML
    """
    try:
        os.mkdir(save_path+'/videos')
    except FileExistsError:
        pass
    
    resultsCDF = Dataset(pathCDF, 'r', format='NETCDF4', parallel=False)

    # Get min and max value for the colorbar
    min_value = np.min(resultsCDF[variable][:].data)
    max_value = np.max(resultsCDF[variable][:].data)

    ## Figure Options ##
    
    fig = plt.figure(figsize=(12,8))
    ax = plt.subplot(111)
    
    ax.set_title(var2str(variable), fontsize=25, pad=15)
    
    # Color bar creation (static one)
    cbar = fig.colorbar(ScalarMappable(cmap=cmap, norm=Normalize(vmin=min_value, vmax=max_value)), ax=ax)
    cbar_ticks = cbar.get_ticks()
    cbar.set_ticks(cbar_ticks)
    
    # Axis tweaks
    ax.tick_params(labelsize=15, direction='in', length=10, width=1, pad=5, color='white')
    ax.set_ylabel(r'y Axis $(km)$', fontsize=15)
    ax.set_xlabel(r'x Axis $(km)$', fontsize=15)
    
    frames = []
    
    for iteration_nb in range(resultsCDF.dimensions['Nt'].size):
    
        # Creating the figure
        myFig = plt.imshow(resultsCDF[variable][:,:,iteration_nb].T, 
                          origin='lower', 
                          cmap=cmap,
                          vmin=min_value,
                          vmax=max_value)
                          
        # Info on elapsed time
        myTitle = 'Elapsed Time = {hour:2d}'.format(hour=int(resultsCDF['t'][iteration_nb]/3600)) + ' hours'
        subtitle = ax.text(0.1,-0.15,myTitle,
                        size=plt.rcParams["axes.titlesize"],
                        ha="center", transform=ax.transAxes)
                        
        # Add the pair to the list of frames
        frames.append([myFig, subtitle])

    resultsCDF.close()
    anim = ArtistAnimation(fig, frames, interval=200, blit=False, repeat=False)
    
    anim.save(save_path+'/videos/'+variable+'.mp4')
    plt.close()
    
    mp4 = open(save_path+'/videos/'+variable+'.mp4','rb').read()
    data_url = "data:"+save_path+"/videos/"+variable+"/mp4;base64," + b64encode(mp4).decode()
    return HTML("""<video width=800 controls>
                         <source src="%s" type="video/mp4">
                   </video>""" % data_url)
                   
def var2str(var_name):
    """Convert names used for computation into better suitables names for plots (specific to the tropopause problem and the examples already implemented)

    :param var_name: name of the variable during computation
    :type var_name: str
    :return: the name of the variable for the plots
    :rtype: str
    """
    switcher_name = {
        'ut':r'$u_t \ (m.s^{-1})$',
        'vt':r'$v_t \ (m.s^{-1})$',
        'theta_t':r"$\theta^{'}_{tp} \ (K)$",
        'us':r'$u_s \ (m.s^{-1})$',
        'vs':r'$v_s \ (m.s^{-1})$',
        'w':r'$w \ (m.s^{-1})$',
        'Delta_z':r'$\Delta z \ (m)$',
        'Delta_T_bb':r'$\Delta T_{bb} \ (K)$',
        'Delta_T_hist':r'$\Delta T_{bb}^{hist} \ (K)$'
    }

    return switcher_name.get(var_name, var_name)

