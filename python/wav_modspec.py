import numpy as np
import scipy.signal

def conv_fft(a, b, shape='full'):

	a = a.ravel()
	sza = len(a)
	
	N = sza + b.shape[0] - 1

	ffta = np.fft.fft(a, N)

	try:
		K = b.shape[1]
	except IndexError:
		K = 1
	ctmp = np.zeros([N, K], dtype=complex)

	if K==1:
		ctmp = np.fft.ifft(ffta*np.fft.fft(b,N), N)
	else:
		for ik in range(K):
			ctmp[:, ik] = np.fft.ifft(ffta*np.fft.fft(b[:, ik], N), N)

	if shape is 'full':
		return ctmp
	elif shape is 'same':
		tocut = (N-sza)/2
		if K==1:
			return ctmp[(1+int(np.ceil(tocut))):(-(int(np.floor(tocut)-1.0)))]
		else:
			return ctmp[ (1+int(np.ceil(tocut))):(-(int(np.floor(tocut)-1.0))), : ]

	return

def cmorlet_wavelet(x, fs, freqs, n=6, normalization=True):

	n_samples = x.shape[0]

	try:
		n_channels = x.shape[1]
	except IndexError:
		n_channels = 1

	n_freqs = len(freqs)

	n_samples_wav = np.round( (2*n/freqs[0])*fs )

	half = np.floor(n_samples_wav/2)

	if (n_samples_wav % 2) == 1:
		time = np.arange(-half, half+1)/fs
	else:
		time = np.arange(-(half-1), half+1)/fs

	wfam = np.zeros([len(time), n_freqs], dtype=complex)

	for iwav in range(n_freqs):
		s = n/(2*np.pi*freqs[iwav])
		gaussian_win = np.exp((-time**2)/(2*s**2))
		sinwave = np.exp(2*np.pi*1j**freqs[iwav]*time)

		if normalization:
			A = 1. / ((s**2) * np.pi) ** (1./4)
		else:
			A = 1.

		wfam[:, iwav] = A * sinwave * gaussian_win

		wcoef = np.zeros([n_samples, n_freqs, n_channels], dtype=complex)

	if n_channels == 1:
		tmp = conv_fft(x, wfam, 'same')
		wcoef[:, :, 0] = tmp		
	else:
		for i_channel in range(n_channels):
			x_tmp = x[:, i_channel]
			tmp = conv_fft(x_tmp, wfam, 'same')
			wcoef[:, :, i_channel] = tmp			

	return wcoef, wfam


def rfft_transform(x, n=None, dim=None):

	shape_x = x.shape
	dim_x = len(shape_x)

	assert dim_x<=2

	dim_def = dim_x-1

	if dim is None:
		dim = dim_def

	if n is None:
		n = shape_x[dim]

	y = np.fft.fft(x, axis=dim)

	if n%2 == 0:     
		n_new = int((n / 2) + 1)
	else:
		n_new = int((n + 1) / 2)

	if dim == 0:
		y_real = y[:, 0:n_new]		
	else:
		y_real = y[0:n_new, :]

	return y_real

def rfft_psd(x, fs, n_fft=None, win_funct = 'blackmanharris'):

	n_samples, n_channels = x.shape

	if n_fft is None:
		n_fft = x.shape[0]

	win = scipy.signal.get_window(win_funct, n_samples)

	win.shape = (n_samples, 1)

	win_rms = np.sqrt(np.sum(np.square(win)) / n_samples)
	win_mat = np.tile(win, n_channels)
	# check if it is working properly

	x = np.multiply(x, win_mat)

	Xt = rfft_transform(x, n_fft)

	Xt = Xt / win_rms

	X_pwr = np.multiply(Xt, np.conj(Xt))
	X_pwr = X_pwr * (1/np.square(n_fft))

	if n_fft % 2 != 0:
		n_freqs = (n_fft + 1) / 2
		X_pwr[1:-1, :] = X_pwr[1:-1, :] * 2
	
	else:
		n_freqs = (n_fft / 2) + 1
		X_pwr[1:-1, :] = X_pwr[1:-1, :] * 2

	rfft_f_delta = (fs / n_fft)

	rfft_psd_res = np.divide(X_pwr, rfft_f_delta)

	n_freqs = int(n_freqs)

	rfft_f_axis = np.asarray(range(0, n_freqs)) * rfft_f_delta			# check

	rfft_result = Xt

	return rfft_result, rfft_f_delta, rfft_psd_res, rfft_f_axis

def wavelet_spectrogram(x, fs):

	freq_vct = np.array(range(1, int(np.floor(fs / 2) + 1)))

	# Time delta
	t_delta = 1 / fs
	
	# Frequency delta
	fs_delta = freq_vct[1] - freq_vct[0]

	# Create time vector 'time_vct' for signal 'x'
	time_vct = np.array(range(0, np.size(x, 0))) / fs

	# Number of samples
	n_samples  = np.size(x, 0)

	# Wavelet transform
	wavelet_coefficients, _ = cmorlet_wavelet(x, fs, freq_vct)

	# Power from Wavelet coefficients
	power_spectrogram = np.square(np.abs(wavelet_coefficients))
	power_spectrogram = power_spectrogram * 2 / (fs * n_samples)

	return wavelet_coefficients, power_spectrogram

def wavelet_modulation_spectrogram(x, fs):

	wavelet_coefficients, power_spectrogram = wavelet_spectrogram(x, fs)
	n_windows, n_freqs, n_channels =  wavelet_coefficients.shape

	fft_factor_x = 1
	n_fft_x =  fft_factor_x * n_windows	

	for i_channel in range(0, n_channels):

		spectrogram_1ch = np.sqrt(power_spectrogram[:, :, i_channel])
    
		rfft_result, rfft_freq_delta, rfft_psd_res, rfft_freq_axis = rfft_psd(spectrogram_1ch, fs, n_fft_x)

		fmod_ax = rfft_freq_axis
    
		fmod_delta = rfft_freq_delta

		n_freqsmod = np.size(fmod_ax)

		if i_channel == 0:
			rfft_modspec = np.zeros([n_freqs, n_freqsmod, n_channels], dtype = complex)
			pwr_modspec  = np.zeros([n_freqs, n_freqsmod, n_channels], dtype = complex)
    
		rfft_modspec[:, :, i_channel] = np.transpose(rfft_result)			# TIREI TRANSPOSE
		pwr_modspec[:, :, i_channel] = np.transpose(rfft_psd_res)			# TIREI TRANSPOSE

	pwr_modspec = pwr_modspec / fmod_delta

	return pwr_modspec, rfft_modspec

if __name__ == '__main__':
	
	fs = 256
	t_5s = np.arange(5*fs)/fs
	freqs = np.arange(1,101)
	x = np.asarray([np.sin(8*2*np.pi*t_5s), np.sin(25*2*np.pi*t_5s)])
	#x = x.ravel()	
	x = np.transpose(x)
	t = np.arange(len(x))/fs

	a = wavelet_modulation_spectrogram(x, fs)

	
