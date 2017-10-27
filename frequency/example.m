clear, clc, close all

% signal parameters
fs = 512;
f0 = 137.365;
T = 0.5;
N = round(T*fs);
SNR = 60;

% signal generation
t = (0:N-1)/fs;
s = sin(2*pi*f0*t);
n = 1/sqrt(2*(10^(SNR/10)))*randn(1, N);
x = s + n;

% signal frequncy estimation
f0hat = frequency(x, fs);

% error of the estimation
mse = (f0hat - f0)/f0*100;

% display the results
f0hatstr = num2str(f0hat);
msestr = num2str(mse);
disp(['Estimated frequency = ' f0hatstr])
disp(['Error of the estimation = ' msestr, ' %'])

commandwindow