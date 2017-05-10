# Amplitude Modulation Toolbox

The Amplitude Modulation Toolbox, for MATLAB, provides functions to compute and visualize the frequency-frequency-domain representation of real signals.

The Toolbox includes a GUI implementation, which facilitates the Amplitude Modulation analysis by allowing changing parameters on line.

In summary, the frequency-frequency representation of a signal is computed by performing two transformations.

<p align="center">
<b>Real signal (time domain)</b>  
</p>
<p align="center">
<i>First Transformation</i>
</p>
<p align="center">
<b>Spectrogram (time-frequency domain)</b>
</p>
<p align="center">
<i>Second Transformation</i>
</p>
<p align="center">
<b>Modulation spectrogram (frequency-frequency domain)</b>
</p>

This Toolbox provides two implementations for the time to time-frequency transformation (*First Transformation*), one based on the STFFT, and the other on wavelet transform (using Complex Morlet wavelet). The time-frequency to frequency-frequency transformation (*Second Transformation*) is carried out with the FFT.

## Examples
Besides the functions to compute and visualize the frequency-frequency representation of real signals, example data and scripts are provided.

### Example 1: `example_01.m`
This example shows the usage and the commands accepted by the high level GUI to explore Amplitude Modulation in ECG and EEG data. The GUI can be called with the functions:
`explore_strfft_am_gui()` which uses STFFT, and `explore_wavelet_am_gui()` based on wavelet transformation. Further details in their use refer to the comments in `example_01.m`.  

![stfft](https://cloud.githubusercontent.com/assets/8238803/25900142/67a297da-3560-11e7-8112-16a7f6c3e637.png)
STFFT-based Amplitude Modulation analysis GUI  
</br>

![wavelet](https://cloud.githubusercontent.com/assets/8238803/25900150/6bf2b93c-3560-11e7-8dd4-084b23c925b5.png)
STFFT-based Amplitude Modulation analysis GUI

### Example 2: `example_02.m`
This script shows presents the details in the usage of the low level functions to carry on the signal transformations, as well as plotting functions. Refer to the comments in `example_02.m`
