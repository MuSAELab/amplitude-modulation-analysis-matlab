"""
Amplitude Modulation Analysis Toolbox
"""

import numpy as np
import scipy.signal

def conv_fft(a, b, shape='full'):
    """
    C = CONV_FFT(A, B) convolves vectors A and a vector or matrix B;
    if B is a vector, C is the convolution between A and B
    if B is a matrix (B is a collection of K kernels), [samples,K]
    C is matrix that  contains the result of the convolution of A
    and each kernel in B. C = [convolution_length, K]
 
    C has the same class as A and B
 
    This function is based in the Convolution Theorem.
    C = ifft(fft(A).*fft(B(:,k)) for k = 1, ... K kernels
 
    The resulting convolution vector is
    length MAX([LENGTH(A)+LENGTH(B)-1,LENGTH(A),LENGTH(B)]).
    If A and B are vectors of polynomial coefficients, convolving
    them is equivalent to  multiplying the two polynomials.

    C = CONV_FFT(A, B, SHAPE) returns a subsection of the convolution with size
    specified by SHAPE:
      'full'  - (default) returns the full convolution,
      'same'  - returns the central part of the convolution
               that is the same size as A.
    """
    # input vector 'a' to 1 dimension
    a = a.ravel()
    sza = len(a)
    
    N = sza + b.shape[0] - 1

    ffta = np.fft.fft(a, N)
    fftb = np.fft.fft(b, N, axis=0)

    # full convolution is computed
    try:
        K = b.shape[1]
    except IndexError:
        K = 1
    ctmp = np.zeros([N, K], dtype=complex)

    if K==1:
        ctmp = np.fft.ifft(ffta * fftb, N)
    else:
        ctmp = np.fft.ifft(ffta[:, np.newaxis] * fftb, N, axis=0 )

    if shape is 'full':
        return ctmp
    elif shape is 'same':
        tocut = (N-sza)/2
        if K==1:
            return ctmp[(int(np.ceil(tocut))):(-(int(np.floor(tocut))))]
        else:
            return ctmp[(int(np.ceil(tocut))):(-(int(np.floor(tocut)))), : ]

    return

def epoching(data, samples_epoch, samples_overlap = 0):
#CHECKED
    """
    Given a 2D array of the shape [n_samples, n_channels]    
    Creates a 3D array of the shape [wlength_samples, n_channels, n_epochs]
    
    Arguments
    data:  2D Data [n_samples, n_channels]
    samples_epoch: Window length in samples
    samples_overlap: Overlap between windows in samples, if it is not specified
                     samples_overlap = 0
    """ 
    try:
        data.shape[1]
    except IndexError:
        data = data[:, np.newaxis]
    
    n_samples = data.shape[0]
    n_channels = data.shape[1]

    samples_shift = samples_epoch - samples_overlap

    n_epochs =  int(np.floor( (n_samples - samples_epoch) / float(samples_shift) ) + 1 )

    #markers indicates where the epoch starts, and the epoch contains samples_epoch rows
    markers = np.asarray(range(0,n_epochs)) * samples_shift;
    markers = markers.astype(int)
    print(markers)
    #Divide data in epochs
    epochs = np.zeros((samples_epoch, n_channels, n_epochs));

    for i_epoch in range(0,n_epochs):
        epochs[:,:,i_epoch] = data[ markers[i_epoch] : markers[i_epoch] + samples_epoch ,:]
        
    print([markers[-1] + samples_epoch])
    if ( (markers[-1] + samples_epoch) < n_samples): 
        remainder = data[markers[-1] + samples_epoch : n_samples, :]
    else:
        remainder = np.asarray([])
    
    return epochs, remainder


def cmorlet_wavelet(x, fs, freqs, n=6, normalization=True):
#CHECKED
    '''
    [wcoef,wfam] = cmorlet_wavelet(x, fs, freqs, n, normalization)
    CMORLET_WAVELET Wavelet transform of ONE signal
       Continuous wavelet tranform using the Complex Morlet Wavelet
    
       WCOEF has the same class as X
     
       WCOEF = CMORLET_WAVELET(SIGNAL, FS, FREQS, N, NORMALIZATION) 
       obtains the Wavelet transform using the Complex Morlet wavelet 
       as mother wavelet using the definition 
       cmor = A * exp(1i*2*pi*f*time)* exp(-time.^2./(2*s^2))
       with s = n / (2*pi*f) 
       
       if NORMALIZATION  is absent or equal to TRUE:
       A = 1 / ( sqrt(szx)*((s(iwav)^2)*pi)^(1/4) );
     
       X     = real-valued signal or set of signals [n_samples, n_channels]
       FS    = Sampling frequency
       FREQS  = Array with the frequencies where the wavelet convolution is computed
       n     = Relation between the number of cicles allowed by the Gaussian
               defaul n = 6 
       NORMALIZATION = Defines the value of A
                       default normalization = true
       WCOEF =  complex coefficients of the wavelet transform for frequencies in FREQS 
                WCOEF = [n_samples, numel(FREQS),n_channels]
       WFAM  =  wavelet family (set of wavlets) used. Size [n_samples, numel(FREQS)]
    '''
    # number of samples
    n_samples = x.shape[0]
    
    # number of channels
    try:
        n_channels = x.shape[1]
    except IndexError:
        n_channels = 1

    # number of wavelets
    n_freqs = len(freqs)

    # number of samples for Wavetet family
    # This is equal to the number of samples needed to represent 2*n cycles 
    # of a sine with frequency = fres(1)[Hz], sampled at fs [Hz]. 
    # This is done to ensure that every wavelet in the wavalet family will be 
    # close to 0 in the negative and positive edges
    n_samples_wav = np.round( (2*n/freqs[0])*fs )

    # create time vector for Wavelet family
    half = np.floor(n_samples_wav/2)
    if (n_samples_wav % 2) == 1:
        time = np.arange(-half, half+1)/fs
    else:
        time = np.arange(-(half-1), half+1)/fs

    # initialize Wavelet family matrix
    wfam = np.zeros([len(time), n_freqs], dtype=complex)

    # for each frequency defined in FREQ, create its respective Wavelet
    for iwav in range(n_freqs):
        s = n/(2*np.pi*freqs[iwav])
        gaussian_win = np.exp((-time**2)/(2*s**2))
        sinwave = np.exp(2*np.pi*1j*freqs[iwav]*time)
        if normalization:
            # each wavelet has unit energy sum(abs(wavelet).^2)) = 1
            A = 1. / ((s**2) * np.pi) ** (1./4)
        else:
            A = 1.
        # Complex Morlet wavelet
        wfam[:, iwav] = A * sinwave * gaussian_win

    wcoef = np.zeros([n_samples, n_freqs, n_channels], dtype=complex)

    if n_channels == 1:
        # one channel
        tmp = conv_fft(x, wfam, 'same')
        wcoef[:, :, 0] = tmp        
    else:
        # convolution between signal X and the each Wavelt in the Wavelet family
        for i_channel in range(n_channels):
            x_tmp = x[:, i_channel]
            tmp = conv_fft(x_tmp, wfam, 'same')
            wcoef[:, :, i_channel] = tmp            

    return wcoef, wfam


def rfft_transform(x, n=None, dim=0):
# CHECKED    
    '''
    Y = RFFT_TRANSFORM(X, n, dim)
    
    RFFT Real Fast Fourier Transform
       The RFFT function returns the DFT of a real signal (or signals) in X,
       using the FFT algorithm.
       The output Y contains ONLY the positive frequencies
    
       The usage for RFFT() is the same as the one for FFT()
       Refer to FFT() help for further details
    
       Considering a real signal A with B = fft(A), B is Hermitian symmetric,
       i.e. B(-1) = conj(B(1)), therefore the complete spectrum B
       can be found by using with only the positive frequencies in B
    
       The RFFT uses the Hermitian Symmetry property, returning ONLY positive frequencies
       In this way the spectrum size has approx. 50% the size of the complete
       spectrum, improving storage without losing information
    
     Example:
    
     n = size(X ,1); %Assuming one or more column-vector signals
     Y = RFFT(X);
     X_recover = IRFFT(Y, n);
    '''
    # shape of x
    shape_x = x.shape
    # number of dimentions
    dim_x = len(shape_x)
    
    # limits to 2-dimention data
    assert dim_x<=2
    
    # check shape of X, and set n and dim defaults
    if dim_x == 1:
        dim_def = 0
    else:
        if shape_x[0] == 1:
            # shape [1, n_samples] (row vector)
            dim_def = 1
        elif shape_x[1] == 1:
            # shape [n_samples, 1] (column vector)
            dim_def = 0 
        else:
            # X is a 2D Matrix, a shape [n_samples, n_channels] is asummed
            dim_def = 0;
    
    if dim is None:
        dim = dim_def
    
    if n is None:
        n = shape_x[dim]
    
    # FFT
    y = np.fft.fft(x, n=n, axis=dim)
    
    # points to keep
    if n%2 == 0:    
        # even case
        n_new = int((n / 2) + 1)
    else:
        # odd case
        n_new = int((n + 1) / 2)
    
    if dim_x == 1:
        y_real = y[0:n_new]
    else:
        if dim == 0:
            y_real = y[0:n_new,:]        
        else:
            y_real = y[:, 0:n_new]
    
    return y_real


def rfft_psd(x, fs, n_fft=None, win_funct = 'blackmanharris', channel_names=None):
# CHECKED
#TODO Add the funtionality of handling the channel_names
    '''
    psd_dict = RFFT_PSD(x, fs, n_fft, win_funct, channel_names)
    
     This function computes the PSD for one or a set of REAL signals 'x'
     rFFT and PSD fields have the same numeric class as 'x'
    
     INPUTS:
     x             Real-valued signal or set of signals [n_samples, n_channels]
     fs            Sampling frequency (Hz)
     Optional:
     n_fft         Number of elements to perform FFT (default: n_samples)
     win_funct     Window to apply to the data in 'x' (default: 'blackmanharris')
     channel_names Cell array with names of channels (default: {'CH1', 'CH2', ... , 'CHn_channels'})
    
     OUTPUTS:
     psd_dict.    Output dictionary
       rFFT          First half of the FFT(x) (u), scaled by the Window RMS
       PSD           Power Spectrum Density (u^2 / Hz)
       fs            Sampling frequency (Hz)
       freq_axis     Frequency axis for rFFT and PSD (Hz)
       freq_delta    Frequency axis step (Hz)
       n_samples     Number of samples of the signal or signals 'x'
       n_fft         Number of elements utilized to perform FFT
       win_funct     Window applied to the data in 'x'
       channel_names Cell array with names of channels
    '''
    
    n_samples, n_channels = x.shape
    
    # validate 'n_fft' argument
    if n_fft is None:
        n_fft = n_samples

    # windowing data
    win = scipy.signal.get_window(win_funct, n_samples, fftbins=False)
    win.shape = (n_samples, 1)
    win_rms = np.sqrt(np.sum(np.square(win)) / n_samples)
    win_mat = np.tile(win, n_channels)
    x = np.multiply(x, win_mat)

    # real FFT with zero padding if n_fft ~= n_samples
    Xt = rfft_transform(x, n_fft)
    # spectrum scaled by window RMS
    Xt = Xt / win_rms
    # power spectrum
    X_pwr = np.multiply(Xt, np.conj(Xt))
    X_pwr = X_pwr * (1/np.square(n_fft))

    # adjust for even and odd number of elements
    if n_fft % 2 != 0:
       # odd case
        n_freqs = (n_fft + 1) / 2
       # double all frequency components except DC component 
        X_pwr[1:, :] = X_pwr[1:, :] * 2
    
    else:
       # even case 
        n_freqs = (n_fft / 2) + 1
       # double all frequency components except DC and fs/2 components
        X_pwr[1:-1, :] = X_pwr[1:-1, :] * 2
    
    # frequency axis step
    f_delta = (fs / n_fft)
    # scale PSD with the frequency step
    psd = np.divide(X_pwr, f_delta)

    # frequency axis for spectrum
    n_freqs = int(n_freqs)
    f_axis = np.asarray(range(0, n_freqs)) * f_delta
    
    # output 'psd_dict' dictionary
    psd_dict = {}
    psd_dict['rFFT'] = Xt
    psd_dict['PSD'] = psd
    psd_dict['fs'] = fs
    psd_dict['freq_axis'] = f_axis
    psd_dict['freq_delta'] = f_delta
    psd_dict['n_samples'] = n_samples
    psd_dict['n_fft'] = n_fft
    psd_dict['win_funct'] = win_funct
    psd_dict['channel_names'] = channel_names
    
    return psd_dict

def strfft_spectrogram():
    pass

def wavelet_spectrogram(x, fs, n_cycles=6, freq_vct=None, channel_names=None):
#CHECKED
#TODO Add the funtionality of handling the channel_names
    '''
    spectrogram_dict = WAVELET_SPECTROGRAM(x, fs, n_cycles, freq_vct, channel_names)
     This function computes the Spectrogram using the Complex Morlet wavelet
     for one or a set of REAL signals 'x'
     rFFT_spectrogram and pwr_spectrogram have the same numeric class as 'x'
    
     INPUTS:
      x             Real-valued signal or set of signals [n_samples, n_channels]
      fs            Sampling frequency (Hz)
     Optional:
      freq_vct      Indicates the frequencies for the Wavelet transform
                     Default = 1 : floor(fs/2);
      n_cycles      Approx. number of cycles in the Gaussian kernel (default = 6)
     channel_names  Cell array with names of channels (default: {'CH1', 'CH2', ... , 'CHn_channels'})
    
     OUTPUTS:
      spectrogram_dict. Output dictionary
       wavelet_coef      Coefficients of the Wavelet transformation (u)
       pwr_spectrogram   Scaled power
       fs                Sampling frequency (Hz);
       freq_axis         Frequency axis for rFFT_ and pwr_spectrogram (Hz)
       freq_delta        Frequency axis step (Hz)
       time_axis         Time axis for rFFT_ and pwr_spectrogram (s)
       time_delta        Time axis step (s)
       n_cycles          Number of cycles in the Gaussion kernel
       wavelet_kernels   Wavelet kernels to obtain the wavelet coefficients
       n_samples         Number of samples of the signal or signals 'x'
       channel_names     Cell array with names of channels
    '''
    # validate 'freq_vct' argument
    if freq_vct is None:
        freq_vct = np.array(range(1, int(np.floor(fs / 2) + 1)))
    
    # Time delta
    t_delta = 1 / fs
    
    # Frequency delta
    f_delta = freq_vct[1] - freq_vct[0]

    # Create time vector 'time_vct' for signal 'x'
    time_vct = np.array(range(0, np.size(x, 0))) / fs

    # Number of samples
    n_samples  = np.size(x, 0)

    # Wavelet transform
    wavelet_coefficients, wavelet_family = cmorlet_wavelet(x, fs, freq_vct)

    # Power from Wavelet coefficients
    power_spectrogram = np.square(np.abs(wavelet_coefficients))
    power_spectrogram = power_spectrogram * 2 / (fs * n_samples)

    # output 'spectrogram_dict' dictionary
    spectrogram_dict = {}
    spectrogram_dict['wavelet_coef'] = wavelet_coefficients
    spectrogram_dict['pwr_spectrogram'] = power_spectrogram
    spectrogram_dict['fs'] = fs
    spectrogram_dict['freq_axis'] = freq_vct
    spectrogram_dict['freq_delta'] = f_delta
    spectrogram_dict['time_axis'] = time_vct
    spectrogram_dict['time_delta'] = t_delta
    spectrogram_dict['n_cycles'] = n_cycles
    spectrogram_dict['wavelet_kernels'] = wavelet_family    
    spectrogram_dict['n_samples'] = n_samples
    spectrogram_dict['channel_names'] = channel_names

    return spectrogram_dict

def strfft_modulation_spectrogram():
    pass


def wavelet_modulation_spectrogram(x, fs, n_cycles=6, freq_vct=None, fft_factor_x=1, win_funct_x='blackmanharris', channel_names=None):

    spectrogram_data = wavelet_spectrogram(x, fs, n_cycles, freq_vct, channel_names)
    n_windows, n_freqs, n_channels =  spectrogram_data['wavelet_coef'].shape

    # number of elements for FFT of the spectrogram
    n_fft_x =  fft_factor_x * n_windows    

    fs_mod = fs

    # the AM analysis is made in the Amplitude derived from the Power Spectrogram
    for i_channel in range(0, n_channels):
        # data to generate the Modulation Spectrogram
        spectrogram_1ch = np.sqrt(spectrogram_data['pwr_spectrogram'][:, :, i_channel])
        # Compute rfft_psd on each frequency timeseries
        psd_dict = rfft_psd(spectrogram_1ch, fs, n_fft_x)
    
        rfft_result = psd_dict['rFFT']
        rfft_psd_res = psd_dict['PSD']
       
        if i_channel == 0:
            # modulation frequency axis
            fmod_ax = psd_dict['freq_axis']
            # modulation frequency delta
            fmod_delta = psd_dict['freq_delta']
            n_freqsmod = np.size(fmod_ax)
            # initialize 'rFFT_modspec'  and 'pwr_modspec'
            rfft_modspec = np.zeros([n_freqs, n_freqsmod, n_channels], dtype = complex)
            pwr_modspec  = np.zeros([n_freqs, n_freqsmod, n_channels], dtype = complex)
    
        # rFFT data
        rfft_modspec[:, :, i_channel] = np.transpose(rfft_result)            
        # power data
        pwr_modspec[:, :, i_channel] = np.transpose(rfft_psd_res)
    
    # scale 'pwr_modspec' by modulation-frequency delta    
    pwr_modspec = pwr_modspec / fmod_delta

    # output 'modspectrogram_dict' dictionary
    modspectrogram_dict = {}
    modspectrogram_dict['rFFT_modspec'] = rfft_modspec
    modspectrogram_dict['pwr_modspec'] = pwr_modspec
    modspectrogram_dict['fs'] = fs
    modspectrogram_dict['fs_mod'] = fs_mod
    modspectrogram_dict['freq_axis'] = spectrogram_data['freq_axis']
    modspectrogram_dict['freq_delta'] = spectrogram_data['freq_delta']
    modspectrogram_dict['modfreq_axis'] = fmod_ax
    modspectrogram_dict['modfreq_delta'] = fmod_delta
    modspectrogram_dict['n_fft_x'] = n_fft_x
    modspectrogram_dict['win_funct_x'] = win_funct_x
    modspectrogram_dict['n_samples'] = spectrogram_data['n_samples'] 
    modspectrogram_dict['channel_names'] = channel_names
    modspectrogram_dict['spectrogram_dict'] = spectrogram_data
    
    return modspectrogram_dict

def plot_modulation_spectrogram_data(modalation_spectrogram_data, ix=None, f_range=None, fmod_range=None, c_range=None):
    pass

def plot_spectrogram_data(spectrogram_data, ix=None, t_range=None, f_range=None, c_range=None):
    pass

def plot_psd_data(psd_data, ix=None, p_range=None):
    pass

def plot_signal(x, fs, name=None):
    """Return a foobang

    Optional plotz says to frobnicate the bizbaz first.
    """
    pass







if __name__ == '__main__':
    
    fs = 256
    t_5s = np.arange(5*fs)/fs
    freqs = np.arange(1,101)
    x = np.asarray([np.sin(8*2*np.pi*t_5s), np.sin(25*2*np.pi*t_5s)])
    #x = x.ravel()    
    x = np.transpose(x)
 
    r = wavelet_modulation_spectrogram(x, fs)
    print('Hi')