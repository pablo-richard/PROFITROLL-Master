import numpy as np

def geostwind(a, b, thetatp, params, z=0, fourier=False, verbose=0):
    
    f       = 1e-4
    theta00 = params['theta_00']
    g       = params['g']
    Ns      = params['N_s']
    Nt      = params['N_t']
    N       = Ns if z>0 else Nt
    pate    = g*(Ns-Nt)/(theta00*Ns*Nt)
    
    Pa, Pb = thetatp.shape
    
    thetatphat = np.fft.fft2(thetatp)
    freqx = np.fft.fftfreq(Pa, a/Pa)
    freqy = np.fft.fftfreq(Pb, b/Pb)
    
    vecFreqX = 2*np.pi*freqx
    vecFreqY = 2*np.pi*freqy
    
    KmatX, KmatY = np.meshgrid(vecFreqX, vecFreqY)
    Kmat = np.sqrt(KmatX**2 + KmatY**2).T
    Kmat[np.where(Kmat==0)]=float('Inf')
    
    Mat = pate/Kmat if (z==0) else pate/Kmat * np.exp(-N*Kmat/f*np.abs(z))
    psihat = thetatphat * Mat
    
    if fourier:
        ug = -np.fft.ifft2(psihat*1j*KmatY.T).real
        vg = np.fft.ifft2(psihat*1j*KmatX.T).real
    
    else:
        psi = np.fft.ifft2(psihat).real
        ug = -(np.roll(psi,-1,1)-np.roll(psi,1,1))/(2*a/Pa)
        vg = (np.roll(psi,-1,0)-np.roll(psi,1,0))/(2*b/Pb)

    return ug, vg

def vertwind(a, b, thetatp, thetatpprev, dt, params, z=0, verbose=0):
    
    f       = 1e-4
    theta00 = params['theta_00']
    g       = params['g']
    Ns      = params['N_s']
    Nt      = params['N_t']
    N       = Ns if z>0 else Nt
    pate    = g*(Ns-Nt)/(theta00*Ns*Nt)
    
    Pa, Pb = thetatp.shape
    
    thetatphat = np.fft.fft2(thetatp)
    thetatpprevhat = np.fft.fft2(thetatpprev)
        
    freqx = np.fft.fftfreq(Pa, a/Pa)
    freqy = np.fft.fftfreq(Pb, b/Pb)
    
    vecFreqX = 2*np.pi*freqx
    vecFreqY = 2*np.pi*freqy
    
    KmatX, KmatY = np.meshgrid(vecFreqX, vecFreqY)
    Kmat = np.sqrt(KmatX**2 + KmatY**2).T
    Kmat[np.where(Kmat==0)]=float('Inf') # set to inf. null wavenumbers
    
    Mat = pate/Kmat if (z==0) else pate/Kmat * np.exp(-N*Kmat/f*np.abs(z))
    psihat     = thetatphat * Mat
    psiprevhat = thetatpprevhat * Mat
    
    Kmat[np.where(np.isinf(Kmat))]=0 # set inf. elements back to 0
    
    thetazhat     = theta00/g*(-np.sign(z))*N*Kmat*psihat
    thetazprevhat = theta00/g*(-np.sign(z))*N*Kmat*psiprevhat

    ug = -np.fft.ifft2(psihat*1j*KmatY.T).real
    vg =  np.fft.ifft2(psihat*1j*KmatX.T).real
    
    thetaz = np.fft.ifft2(thetazhat).real
    thetazprev = np.fft.ifft2(thetazprevhat).real
    dtthetaz = (thetaz-thetazprev)/dt
    dxthetaz = np.fft.ifft2(thetazhat*1j*KmatX.T).real
    dythetaz = np.fft.ifft2(thetazhat*1j*KmatY.T).real
    
    w = -dtthetaz - ug*dxthetaz - vg*dythetaz
    w *= g/(N**2 * theta00)
    
    return w
