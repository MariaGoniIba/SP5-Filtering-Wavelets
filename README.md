# Signal Processing (SP) project 5: FIR filtering and Wavelets filtering
Filtering a signal both using FIR filters and Wavelets.

# Credit
The dataset and proposal of the exercise is from the Udemy course [Signal Processing Problems, solved in Matlab and in Python](https://www.udemy.com/course/signal-processing/). I highly recommend this course for those interested in signal processing.

# Proposal
We are given a noisy signal '_signal_' and the time series of the filtered signal '_signalFIR_' using a FIR filter and the filtered signal '_signalMV_' using Morlet Wavelets. 
The goal is to find out how to create a FIR filter and a Morlet wavelet to generate similar results.

# Procedure
First, we plot the given signals both in the time and frequency domain, to be able to set up the cut-off frequencies of the filters.
<p align="center">
    <img width="600" src="https://github.com/MariaGoniIba/SP5-Filtering-Wavelets/blob/main/ResultsToObtain.jpg">
</p>

## FIR filter
I apply a band-pass filter at frequencies between 11 and 15 Hz. I study the frequency response of the filter before applying it.
<p align="center">
    <img width="600" src="https://github.com/MariaGoniIba/SP5-Filtering-Wavelets/blob/main/FilterResponse.jpg">
</p>

## Wavelets filter
We create the Morlet Wavelet.
```
% centered time vector
wavtime = -3:1/srate:3;

% parameters
freq = 13; % peak frequency
fwhm = 0.15; % full-width at half-maximum in seconds
csw  = cos(2*pi*freq*wavtime); % cosine wave
gaussian = exp( -(4*log(2)*wavtime.^2) / fwhm^2 ); % Gaussian
% Morlet wavelet
MorletWavelet = csw .* gaussian;

plot(MorletWavelet)
set(gca,'xlim',[5000 7500])
```
<p align="center">
    <img width="600" src="https://github.com/MariaGoniIba/SP5-Filtering-Wavelets/blob/main/Wavelet.jpg">
</p>

Thinking back to the convolution theorem from a [previous project](https://github.com/MariaGoniIba/SP4-Frequency-domain-mean-smoothing-filter-Convolution) and implementing convolution as spectral multiplication, we can imagine why Morlet Wavelet is useful for narrow band filtering a time domain signal.

If the filtered signal is larger than the original signal, this is a problem of normalization. 
This is because of the gain that you have in the frequency domain, you will be amplyfing it in the signal. 
But normalising the wavelet in the time domain to get the maximum value of 1 in the frequency domain turns out to be very tricky. 
So the best way to do it is to normalize in the frequency domain, dividing the entire spectrum by the maximum value. 
Therefore, we don’t want to use the matlab _conv_ function, but implement convolution manually. 

Remember that the length of the result of the convolution is the number of points in the signal plus the number of points in the kernel -1. 
When you implement the convolution manually the order doesn’t matter because multiplication is cumulative. 
But if you use the matlab conv function it matters, because matlab assumes that the first is signal and the second the kernel.

Therefore:
```
% Manual convolution
nConv = N + length(wavtime) -1;
halfw = floor(length(wavtime)/2)+1;

% spectrum of wavelet
morwavX = fft(MorletWavelet,nConv);

% now normalize in the frequency domain
morwavX = morwavX ./ max(abs(morwavX));

% rest of convolution
convres = ifft( morwavX .* fft(signal',nConv) );
convres = real( convres(halfw:end-halfw+1) );

% power
yfilt2f = abs(fft(convres)).^2;
```

## Solution
We finally plot in the time and frequency-domain the results of both filters.
<p align="center">
    <img width="600" src="https://github.com/MariaGoniIba/SP5-Filtering-Wavelets/blob/main/Results.jpg">
</p>



