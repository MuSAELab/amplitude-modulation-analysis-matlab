# Amplitude Modulation Analysis Toolbox

The Amplitude Modulation Analysis Toolbox, for **MATLAB / Octave**, provides functions to compute and visualize the time-frequency and frequency-frequency domain representations of real signals. The **Python 3** version of this toolbox can be found here: https://github.com/MuSAELab/amplitude-modulation-analysis-module

The Toolbox includes a GUI implementation, which facilitates the exploration of the amplitude modulation by allowing changing parameters online.

In summary, the frequency-frequency representation of a signal is computed by performing two transformations:

1. Time domain (real signal) to time-frequency (spectrogram)
2. Time-frequency (spectrogram) to frequency-frequency (modulation spectrogram)

![diagram](https://user-images.githubusercontent.com/8238803/32670991-b1a4e542-c613-11e7-8408-bdc1cc3e0bf1.png)
Signal processing steps involved in the calculation of the modulation spectrogram from the amplitude spectrogram of signal. The block |abs| indicates the absolute value, and the FT indicates the use of the Fourier transform.

This Toolbox provides two implementations for the time to time-frequency transformation, one based on the short-time FFT (STFFT) and the other on the continuous wavelet transform (CWF) using the Complex Morlet wavelet. The time-frequency to frequency-frequency transformation is carried out with the FFT.

## Examples
Besides the functions to compute and visualize the frequency-frequency representation of real signals, example data and scripts are provided.

### Example 1: `example_01.m`
This example shows the usage and the commands accepted by the high level GUI to explore Amplitude Modulation in ECG and EEG data. The GUI can be called with the functions:
`explore_strfft_am_gui()` which uses STFFT, and `explore_wavelet_am_gui()` based on wavelet transformation. Further details in their use refer to the comments in `example_01.m`.  

![stfft](https://cloud.githubusercontent.com/assets/8238803/25900142/67a297da-3560-11e7-8112-16a7f6c3e637.png)
STFFT-based amplitude modulation analysis GUI  
</br>

![wavelet](https://cloud.githubusercontent.com/assets/8238803/25900150/6bf2b93c-3560-11e7-8dd4-084b23c925b5.png)
CWT-based amplitude modulation analysis GUI

### Example 2: `example_02.m`
This script shows presents the details in the usage of the low level functions to carry on the signal transformations, as well as plotting functions. Refer to the comments in `example_02.m`

### Citation
If you use this software, please cite the companion [conference article](https://doi.org/10.1109/GlobalSIP45357.2019.8969210).

**BibTeX**
```
@inproceedings{2019_cassani_amAn_Opensource_2019,
  author = {Cassani, Raymundo and Albuquerque, Isabela and Monteiro, João and Falk, Tiago H.},
  booktitle = {2019 IEEE Global Conference on Signal and Information Processing (GlobalSIP)},
  doi = {10.1109/GlobalSIP45357.2019.8969210},
  pages = {1--4},
  title = {{AMA: An Open-source Amplitude Modulation Analysis Toolkit for Signal Processing Applications}},
  year = {2019}
}
```
In addition, the conference poster related to the toolbox can be found [here](https://www.castoriscausa.com/files/cassani_2019_amaposter.pdf).


### Acknowledgement
The research is based upon work supported by the Office of the Director of National Intelligence (ODNI), Intelligence Advanced Research Projects Activity (IARPA), via IARPA Contract N°2017 - 17042800005 . The views and conclusions con - tained herein are thos e of the authors and should not be interpreted as necessarily representing the official policies or endorsements, either expressed or implied, of the ODNI, IARPA, or the U.S. Government. The U.S. Government is authorized to reproduce and distribute reprint s for Governmental purposes notwithstanding any copyright annotation thereon.
