%Very crude continuous wavelet transform. A few years ago I tried to understand the
%literature, and this code is the result of that process. There are a lot
%of things I didn't figure out in the end, but the code seems to work.
%Big thanks to Stuart Baker and Stefane Aguiar at Newcastle University for helping me with this!

%This primitive wavelet transform is ideal for low-grade microcontrollers,
%where you have to check for the presence of a known signal.
%It also may be useful for time-frequency analysis, where FFT would fail.

%This does not require any external libraries, and works with both Matlab
%and octave.
clear all;
close all;
clc;

%% general conditions
sampling_rate = 1000; %Hz
signal_length = 2.5; %seconds
%% generate my test signal
sampling_period = 1/sampling_rate;
time = (0:sampling_period:signal_length);
%This is my test signal:
a1 = 1; %relative amplitude of signal 1.
f1 = 25; %relative freqency of signal 1. Hz.
onset1 = 0.3 * sampling_rate; %This is the onset time for the first signal.
a2 = 1; %relative amplitude of signal 2.
f2 = 4; %relative frequency of signal 2. Hz.
onset2 = 1.222 * sampling_rate; %Second signal's onset.
component1 = a1*sin(2*pi*f1*time); %define signals


component1(1:onset1) = 0; %Make sure the signal is only starting at the specified time.
component2 = a2*sin(2*pi*f2*time);


component2(1:onset2) = 0; %Make sure the signal is only starting at the specified time.

%Add some extra chopoffs, totally arbitrarily.
component1(600:2000) = 0;
component1(1800:sampling_rate * signal_length) = 0;

component2(1664:2004) = 0; %Add chopoff in component 2.

test_signal = component1 + component2;


%% Generate scale range and loop the whole lot through it.
%TODO: NEED TO FIGURE OUT HOW THE SCALES RELATE TO FREQUENCY!!!!!!
low_freq_limit = 1; %Hz.
high_freq_limit = 50; %Hz
%resolution = (high_freq_limit - low_freq_limit); %This is 1 Hz per line.
resolution = 500; %how fine each line should be.
base_wavelet_length = 4; %4 seems to work the best. Probably will need to change with the wavelet.
wavelet_natural_frequency = sampling_rate/base_wavelet_length; %As our wavelet is made up from many samples, scale = 1 is this frequency.
highest_wavelet_scale = wavelet_natural_frequency / low_freq_limit; % Low frequency -> high scale
lowest_wavelet_scale = wavelet_natural_frequency / high_freq_limit; % High frequency -> low scale


if(lowest_wavelet_scale < 1)
    error('Highest frequency limit can''t be more than the natural frequency of the wavelet!')
end
scale_range = linspace(lowest_wavelet_scale, highest_wavelet_scale, resolution); %Standard linear scale. Frequency output is NOT linear.


%% Do the actual transform.
for(i = 1:resolution)
    %% create the wavelet.
    wavelet_scale = scale_range(i);  
    
    wavelet_time = linspace(-base_wavelet_length/2, base_wavelet_length/2, wavelet_scale*(base_wavelet_length*2));
    %shannon_wavelet = sinc(wavelet_time ./ 2) .* cos((3*pi*wavelet_time) ./ 2);
    wavelet_real = sinc(wavelet_time ./ 2) .* cos((3*pi*wavelet_time) ./ 2); %This is the Shannon wavelet, real bit.
    wavelet_imag = sinc(wavelet_time ./ 2) .* cos((3*pi*wavelet_time) ./ 2 + pi/2); %This is the Shannon wavelet, imagniary bit, effectively phase-shifted.
    
    %Amplify the wavelet with 1/(scale).
    %Basically, the longer the wavelet, the smaller it needs to be to
    %maintain the same energy.
    wavelet_real = wavelet_real .* (1/(wavelet_scale));
    wavelet_imag = wavelet_imag .* (1/(wavelet_scale));


    %% Apply wavelet transform by convolving the scale with the original signal

    raw_output_real = conv(test_signal, wavelet_real, 'same'); %This signal is AC, and the offsets are compensated for.
    %rectified_output = abs(raw_output_real); %This just takes the real numbers into account.
    %Do the same with the complex numbers too.
    raw_output_imag = conv(test_signal, wavelet_imag, 'same'); %This signal is AC, and the offsets are compensated for.
    rectified_output = sqrt(raw_output_real .^2 + raw_output_imag .^2); %This signal has no negative values any more. This calculates the magnitudes.

    rectified_image(i, :) = (rectified_output); %assemble the rectified output as an image.
end
%normalising the image between 0 and 1.
%For visualisation only!
normalised_image = rectified_image / max(max(rectified_image));

%% Plot.
subplot(2, 1, 1)
plot(test_signal);
title('Input signal');
xlabel('Time [samples]')
axis tight;
subplot(2, 1, 2)
imshow(normalised_image, 'DisplayRange', [0, 1]);
title('Normalised scalogram of the signal [0...1]')
colormap(gca, 'jet')
xlabel('Time [samples]')
ylabel(sprintf('Wavelet scale [%d...%d]', highest_wavelet_scale, lowest_wavelet_scale));

% NOTE: The scalogram does not actually have information about the spectral
% components of the signal, rather, it has information about the energy of
% the scaled wavelet. If you chose your wavelet wisely, you can use it for
% time-frequency analyisis. But that's just one of many applications.