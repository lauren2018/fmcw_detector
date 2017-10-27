%% 
set(0,'DefaultFigureWindowStyle','docked')
%% written by Taehwan Kim (taehwan@berkeley.edu)
format long;
clear; close all; clc;

%% open simulink model
open_system('./model/symodel');

%% set constants
c = 3e8;
timestep = 1e-9;
stoptime = 10e-6;
Fs = 1/timestep;
alpha = 1e15;
omega0 = 0;
phi0 = 0;

fnoise = 10e6;
snoise = 0.00;
tauc=1/(pi*fnoise);

tau_mzi = 100e-9;
tau_match = 100e-9;
reff=alpha*tau_mzi;
gain=1e9;
amplitude=1;

%% run simulation
distance = 20;
tau = 2*distance/c;
% distance = 1:3:16;

dist_res = (0.5/stoptime)*3e8/alpha;
reference_spower = stoptime * exp(-2*tau/tauc) / pi;

N = round(stoptime/timestep);
periods = (1./(tau_mzi*alpha))/timestep;
rightlength = floor(periods.*floor(N./periods));

iter = 100;
padm = 0;
x = zeros(iter, (1+2*padm)*N);
x_windowed = zeros(iter, (1+2*padm)*N);
per = zeros(iter, (1+2*padm)*N/2+1);

freq = 0:Fs/N:Fs/2;

for i=1:iter
    
    fseed = i;
    sseed = 2*i+1;
    sim('symodel');

    x(i,:)=transpose(padarray(vout8.Data(1:end-1),padm*N,0,'both'));
    x_windowed(i,:)=transpose(padarray(cat(1,vout8.Data(1:rightlength),zeros(N-rightlength,1)),padm*N,0,'both'));

    tempx = squeeze(x(i,:));
    tempx_windowed = squeeze(x_windowed(i,:));

    [per(i,:),F] = periodogram(tempx_windowed,hann(length(tempx_windowed)),length(tempx_windowed),Fs);

    fprintf('Just finished iteration #%d\n', i);
    
end

per_avg = mean(per,1);

% per_avg = zeros(avgnum,N/(2*avgnum)+1);
%     for k=1:avgnum
%         [per_avg(k,:),F] = periodogram(tempx_windowed(((k-1)*end/avgnum)+1:((k)*end/avgnum)),hann(length(tempx_windowed(1:end/avgnum))),length(tempx_windowed(1:end/avgnum)),Fs);          
%     end
%     per_avg = mean(per_avg,1);
%     per_avg = per_avg';
%     [Y, I] = max(per_avg);
% save('data.mat','x');
    
%% plot the periodogram-estmated psd for 

% 
% figure(2);
% plot(distance, dist_est_periodogram-distance');
% hold on;
% plot(distance, dist_est_periodogram_correct-distance');


% figure(1);
% x=vout.Data(1:end-1);
% N = length(data_1);
% 
% psdx_1 = periodogram(data_1,hann(N),N,Fs);
% plot(freq,10*log10(psdx_1));
% hold on
% psdx_2 = periodogram(data_2,hann(N),N,Fs);
% plot(freq,10*log10(psdx_2));
% 
% title('Photocurrent PSD');
% axis([0 freq(end) -140 -40]);
% xlabel('Frequency (Hz)');
% grid on;
% ylabel('PSD (dB/Hz)');

% figure(2);
% xdft_correct = fft(x_correct(1,:));
% xdft_correct = xdft_correct(1:N/2+1);
% psdx_correct = (1/(Fs*N)) * abs(xdft_correct).^2;
% psdx_correct(2:end-1) = 2*psdx_correct(2:end-1);
% 
% plot(freq,10*log10(psdx_correct));
% 
% title('Photocurrent PSD');
% axis([0 freq(end) -140 -40]);
% xlabel('Frequency (Hz)');
% grid on;
% ylabel('PSD (dB/Hz)');
% 
% plot(freq,S);
% 
% title('MUSIC');
% xlim([0 freq(end)]);
% xlabel('Frequency (Hz)');
% grid on;