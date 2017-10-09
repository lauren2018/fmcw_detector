%% written by Taehwan Kim (taehwan@berkeley.edu)
clear; close all; clc;
format long;
set(0,'DefaultFigureWindowStyle','docked')

%% open simulink model
open_system('./model/sunlight');

%% set constants
q = 1.6e-19;
Res = 1;
Prx = 1e-3;
c = 3e8;
timestep = 1e-9;
stoptime = 10e-6;
Fs = 1/timestep;
chirpbw = 10e9;
alpha = chirpbw/stoptime;
omega0 = 0;
phi0 = 0;
prx=0;
plo=1/sqrt(2);

N = round(stoptime/timestep);
fnoise = 10e6;
% snoise = q/Res/Prx;
snoise = 1e-14;
tauc=1/(pi*fnoise);
amplitude=1;
freq = 0:Fs/N:Fs/2;

%% feed-forward phase noise correction parameters
tau_mzi = 1e-9;
tau_match = 1e-9;
reff=alpha*tau_mzi;
gain=1e9;

%% run simulation
% distance = 20;
distance = 20;
iter = 1;
avgnum = 4;
padm = 0;

% distance = 1:3:21;

dist_res = (0.5/stoptime)*3e8/alpha;
% distance = [3 4];
tau_list = 2*distance/c;
reference_spower = stoptime * exp(-2*tau_list/tauc) / pi;

periods = (1./(tau_list*alpha))/timestep;
rightlength = floor(periods.*floor(N./periods));

x = zeros(length(distance), iter, (1+2*padm)*N);
x_windowed = zeros(length(distance), iter, (1+2*padm)*N);
x_windowed_correct = zeros(length(distance), iter, (1+2*padm)*N);
signal_est_periodogram = zeros(length(distance), iter);
signal_est_periodogram_avg = zeros(length(distance), iter);
dist_est_periodogram = zeros(length(distance), iter);
dist_est_periodogram_avg = zeros(length(distance), iter);
dist_est_periodogram_lse = zeros(length(distance), iter);
dist_est_periodogram_avg_lse = zeros(length(distance), iter);
dist_est_periodogram_mle = zeros(length(distance), iter);
dist_est_periodogram_avg_mle = zeros(length(distance), iter);
dist_est_music = zeros(length(distance), iter);

for j=1:length(distance)
    for i=1:iter
        tau = tau_list(j);
        fseed = i;
        sseed = 2*i+1;
        sim('sunlight');
        
        datatmp=squeeze(vout.Data);
        
        x(j,i,:)=padarray(datatmp(1:end-1),padm*N,0,'both');
        
        x_windowed(j,i,:)=padarray(cat(1,datatmp(1:rightlength(j)),zeros(N-rightlength(j),1)),padm*N,0,'both');
        
        tempx = squeeze(x(j,i,:));
        tempx_windowed = squeeze(x_windowed(j,i,:));

        [per,F] = periodogram(tempx_windowed,hann(length(tempx_windowed)),length(tempx_windowed),Fs);
        [Y, I] = max(per);
        signal_est_periodogram(j,i) = Y;
        dist_est_periodogram(j,i) = F(I)*c/alpha/2;
        
        solution = pnoise_fit(10*log10(per),F,[N/4 tauc*1e9]);
        dist_est_periodogram_lse(j,i) = c*solution(1)*(F(2)-F(1))/alpha/2;
        
        solution = mlelorentz(per,F,[N/4 tauc*1e9]);
        dist_est_periodogram_mle(j,i) = c*solution(1)*(F(2)-F(1))/alpha/2;
        
        per_avg = zeros(avgnum,N/(2*avgnum)+1);
        for k=1:avgnum
            [per_avg(k,:),F] = periodogram(tempx_windowed(((k-1)*end/avgnum)+1:((k)*end/avgnum)),hann(length(tempx_windowed(1:end/avgnum))),length(tempx_windowed(1:end/avgnum)),Fs);          
        end
        per_avg = mean(per_avg,1);
        per_avg = per_avg';
        [Y, I] = max(per_avg);
        signal_est_periodogram_avg(j,i) = Y;
        dist_est_periodogram_avg(j,i) = F(I)*c/alpha/2;
        
        solution = pnoise_fit(10*log10(per_avg),F,[N/4/avgnum tauc*1e9]);
        dist_est_periodogram_avg_lse(j,i) = c*solution(1)*(F(2)-F(1))/alpha/2;
        
        solution = mlelorentz(per_avg,F,[N/4 tauc*1e9]);
        dist_est_periodogram_avg_mle(j,i) = c*solution(1)*(F(2)-F(1))/alpha/2;

        [S, F] = pmusic(tempx,6,freq,Fs);
        [Y, I] = max(S);
        dist_est_music(j,i) = F(I)*c/alpha/2;

        fprintf('Just finished iteration #%d\n', i);
    end
end
% save('data.mat','x');
    
%% plot the periodogram-estmated psd for 

% data_1 = squeeze(x(1,1,:));
% data_2 = squeeze(x(2,1,:));
% 
% figure(1);
% hold on;
% plot(distance,10*log10(signal_est_periodogram(:,1)));
% plot(distance,10*log10(signal_est_periodogram_correct(:,1)));
% plot(distance, 10*log10(reference_spower));
% plot(distance, ones(1,length(distance))*10*log10(2*tauc));

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