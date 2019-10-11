%% FT test

clear all;
clc;

dt = 1/100; % sampling frequency
T = 100;    % total record time
t = [0:dt:T];
N = length(t);

% define input signal
f1 = 3*pi/180; 
f2 = 25*pi/180;

A1 = 0.7; 
A2 = 2;
sin1 = A1*sin(2*pi*f1*t);
sin2 = A2*sin(2*pi*f2*(t - 23*dt));

input = sin1 + sin2;

noise_var = 0;
noise = randn(1,N)*sqrt(noise_var);

signal = input + noise;

%plot input signal

figure()
hold on
plot(t, signal);
plot(t, input);
xlabel('time (in s)');
ylabel('Signal');
legend('measured input signal', 'ideal input signal');

%% FFT

Fs = 1/dt;  % sampling frequency in Hz
FFT = fft(signal);
FFT_spec2 = abs(FFT/N);
FFT_spec1 = FFT_spec2(1:N/2+1);
FFT_spec1(2:end - 1) = 2*FFT_spec1(2:end - 1);

f = Fs*(0:(N/2))/N;

figure()
hold on
plot(f, FFT_spec1)
xlabel('Frequency (Hz)')
ylabel('|FFT|')
title('FFT')



