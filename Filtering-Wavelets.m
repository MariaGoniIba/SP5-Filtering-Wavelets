clear; clc; close all

load('wavelet_codeChallenge.mat')

N = length(signal);
time  = (0:N-1)/srate;

%% Results to obtain
%Time domain
subplot(3,1,1)
plot(time, signal)
title('Original signal')
xlabel('Time (s)')
subplot(3,1,2)
plot(time, signalFIR)
hold on
plot(time, signalMW, 'r')
title('Filtered signals')
xlabel('Time (s)')
legend('FIR', 'MW')

% Frequency domain
spow = abs( fft( signal )/N ).^2;
sfiltpow = abs( fft( signalFIR )/N ).^2;
sMVpow = abs( fft( signalMW )/N ).^2;
hz = linspace(0,srate/2,floor(N/2)+1);

subplot(3,1,3)
plot(hz,sfiltpow(1:length(hz)))
hold on
plot(hz,sMVpow(1:length(hz)),'r')
title('Amplitude spectra')
xlabel('Frequency (Hz)')
legend('FIR', 'MW')
set(gca,'xlim',[0 20])


%% FIR filter

frange = [ 11 15 ];
transw  = .09;
order   = round( 21*srate/frange(1) );
shape   = [ 0 0 1 1 0 0 ];
frex    = [ 0 frange(1)*(1-transw) frange frange(2)+frange(2)*transw srate/2 ] / (srate/2);
  
filtkern1 = firls(order,frex,shape);

% Explore response kernel
% power spectrum of the filter kernel
filtpow1 = abs(fft(filtkern1)).^2;
% compute the frequencies vector and remove negative frequencies
hz1 = linspace(0,srate/2,floor(length(filtkern1)/2)+1);
filtpow1 = filtpow1(1:length(hz1));

figure(2), clf
subplot(121)
plot(filtkern1,'linew',2)
xlabel('Time points')
title('Filter kernel 1 (firls)')
axis square
 
% plot amplitude spectrum of the filter kernel
subplot(122), hold on
plot(hz1,filtpow1,'ks-','linew',2,'markerfacecolor','w')
plot(frex*(srate/2),shape,'ro-','linew',2,'markerfacecolor','w')
 
% make the plot look nicer
set(gca,'xlim',[0 frange(2)*1.5])
xlabel('Frequency (Hz)'), ylabel('Filter gain 1')
legend({'Actual';'Ideal'})
title('Frequency response of filter 1 (firls)')
axis square

% apply filter
yfilt1 = filtfilt(filtkern1,1,signal);
yfilt1f = abs(fft(yfilt1)).^2;


%% Morlet wavelet

% centered time vector
wavtime = -3:1/srate:3;

% parameters
freq = 12.5; % peak frequency
fwhm = 0.15; % full-width at half-maximum in seconds
csw  = cos(2*pi*freq*wavtime); % cosine wave
gaussian = exp( -(4*log(2)*wavtime.^2) / fwhm^2 ); % Gaussian
% Morlet wavelet
MorletWavelet = csw .* gaussian;

%plot(MorletWavelet)
%set(gca,'xlim',[5000 7500])

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


%% Plotting

% Time-domain
figure(3)
subplot(2,1,1)
plot(time, yfilt1)
hold on
plot(time,convres,'r')
title('Filtered signals')
xlabel('Time (s)')
legend('FIR', 'MW')

% Frequency-domain
hz1 = linspace(0,srate/2,floor(length(yfilt1f)/2)+1);
subplot(2,1,2)
plot(hz1,yfilt1f(1:length(hz1)))
hold on
plot(hz1,yfilt2f(1:length(hz1)),'r')
title('Amplitude spectra')
xlabel('Frequency (Hz)')
legend('FIR', 'MW')
set(gca,'xlim',[0 20])



