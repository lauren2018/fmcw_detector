%% written by Taehwan Kim (taehwan@berkeley.edu)
clear; close all; clc;
format long;
set(0,'DefaultFigureWindowStyle','docked')

%% open simulink model
% open_system('./model/symodel');

%% set constants
q = 1.6e-19;
Res = 1;
c = 3e8;
timestep = 0.25e-9;
stoptime = 10e-6;
Fs = 1/timestep;
chirpbw = 15e9;
alpha = chirpbw/stoptime;
omega0 = 0;
phi0 = 0;
prx=1/sqrt(2);

N = round(stoptime/timestep);
freq = 0:Fs/N:Fs/2;
fnoise = 1e6;
tauc=1/(pi*fnoise);

dc = (c/(4*pi*fnoise))*log(stoptime*pi*fnoise);
% Prx = 10.^(-10:0.5:-9);
% Prx = 0.316227766016840e-9;
% Prx = 10^-10;
Prx = 10^-9;
% snoise_list = q/Res./Prx;
snoise_list =0;
% snoise = 0;
amplitude=1;

%% feed-forward phase noise correction parameters
tau_mzi = 1e-9;
tau_match = 1e-9;
reff=alpha*tau_mzi;
gain=1e9;

%% run simulation
% distance = 20;
% distance = 20;
iter = 1;
avgnum = 10;
padm = 0;
% 
% % distance=100;
distance = 1:3:22;
% 
dist_res = 0.5*c/chirpbw;
% % distance = [3 4];
tau_list = 2*distance/c;
reference_spower = stoptime * exp(-2*tau_list/tauc) / pi;
% 
periods = (1./(tau_list*alpha))/timestep;
rightlength = floor(periods.*floor(N./periods));
% 
x = zeros(length(distance), iter, (1+2*padm)*N);
x_correct = zeros(length(distance), iter, (1+2*padm)*N);
x_windowed = zeros(length(distance), iter, (1+2*padm)*N);
x_windowed_correct = zeros(length(distance), iter, (1+2*padm)*N);
signal_est_periodogram = zeros(length(distance), iter);
signal_est_periodogram_avg = zeros(length(distance), iter);
signal_est_periodogram_correct = zeros(length(distance), iter);

dist_est_periodogram = zeros(length(Prx),length(distance), iter);
dist_est_periodogram_avg = zeros(length(Prx),length(distance), iter);
dist_est_periodogram_lse = zeros(length(Prx),length(distance), iter);
dist_est_periodogram_avg_lse = zeros(length(Prx),length(distance), iter);
dist_est_periodogram_mle = zeros(length(Prx),length(distance), iter);
dist_est_periodogram_avg_mle = zeros(length(Prx),length(distance), iter);
dist_est_music = zeros(length(Prx),length(distance), iter);
dist_est_periodogram_correct = zeros(length(Prx),length(distance), iter);
dist_est_music_correct = zeros(length(Prx),length(distance), iter);

std_est_periodogram = zeros(length(Prx),length(distance));
std_est_periodogram_lse = zeros(length(Prx),length(distance));
std_est_periodogram_mle = zeros(length(Prx),length(distance));

for k=1:length(Prx)
    for j=1:length(distance)
        for i=1:iter
            snoise = snoise_list(k);
            tau = tau_list(j);
            fseed = i;
            sseed = 2*i+1;
            sim('symodel');

            x(j,i,:)=padarray(vout.Data(1:end-1),padm*N,0,'both');
            x_correct(j,i,:)=padarray(vout1.Data(1:end-1),padm*N,0,'both');

            x_windowed(j,i,:)=padarray(cat(1,vout.Data(1:rightlength(j)),zeros(N-rightlength(j),1)),padm*N,0,'both');
            x_windowed_correct(j,i,:)=padarray(cat(1,vout1.Data(1:rightlength(j)),zeros(N-rightlength(j),1)),padm*N,0,'both');

            tempx = squeeze(x(j,i,:));
            tempx_correct = squeeze(x_correct(j,i,:));
            tempx_windowed = squeeze(x_windowed(j,i,:));
            tempx_windowed_correct = squeeze(x_windowed_correct(j,i,:));

            [per,F] = periodogram(tempx_windowed,hann(length(tempx_windowed)),length(tempx_windowed),Fs);

            [Y, I] = max(per);
            signal_est_periodogram(j,i) = Y;
            dist_est_periodogram(k,j,i) = F(I)*c/alpha/2;        

            solution_lse = pnoise_fit(per,F,[N/4 tauc snoise]);
            dist_est_periodogram_lse(k,j,i) = c*solution_lse(1)*(F(2)-F(1))/alpha/2;

            solution_mle = mlelorentz(per,F,N/4, [tauc snoise]);
            dist_est_periodogram_mle(k,j,i) = c*solution_mle(1)*(F(2)-F(1))/alpha/2;

            per_avg = zeros(avgnum,N/(2*avgnum)+1);
            for l=1:avgnum
                [per_avg(l,:),F] = periodogram(tempx_windowed(((l-1)*end/avgnum)+1:((l)*end/avgnum)),hann(length(tempx_windowed(1:end/avgnum))),length(tempx_windowed(1:end/avgnum)),Fs);          
            end
            per_avg = mean(per_avg,1);
            per_avg = per_avg';
            [Y, I] = max(per_avg);
            signal_est_periodogram_avg(j,i) = Y;
            dist_est_periodogram_avg(k,j,i) = F(I)*c/alpha/2;

            solution_avg_lse = pnoise_fit(per_avg,F,[N/4/avgnum tauc snoise]);
            dist_est_periodogram_avg_lse(k,j,i) = c*solution_avg_lse(1)*(F(2)-F(1))/alpha/2;

            solution_avg_mle = mlelorentz(per_avg,F,N/4/avgnum, [tauc snoise]);
            dist_est_periodogram_avg_mle(k,j,i) = c*solution_avg_mle(1)*(F(2)-F(1))/alpha/2;

            [S, F] = pmusic(tempx,6,freq,Fs);
            [Y, I] = max(S);
            dist_est_music(k,j,i) = F(I)*c/alpha/2;

            [per_correct,F] = periodogram(tempx_windowed_correct,hann(length(tempx_windowed_correct)),length(tempx_windowed_correct),Fs);
            [Y, I] = max(per_correct);
            signal_est_periodogram_correct(j,i) = Y;
            dist_est_periodogram_correct(k,j,i) = F(I)*c/alpha/2;

            [S_correct, F] = pmusic(tempx_correct,6,freq,Fs);
            [Y, I] = max(S_correct);
            dist_est_music_correct(k,j,i) = F(I)*c/alpha/2;
            fprintf('Just finished iteration #%d-%d\n', k,i);
        end
        std_est_periodogram(k,j) = std(dist_est_periodogram(k,j,:));
        std_est_periodogram_lse(k,j) = std(dist_est_periodogram_lse(k,j,:));
        std_est_periodogram_mle(k,j) = std(dist_est_periodogram_mle(k,j,:));
    end
end
% save('data.mat','x');
    
%% plot the periodogram-estmated psd for 

% data_1 = squeeze(x(1,1,:));
% data_2 = squeeze(x(2,1,:));

figure(1);
hold on;
plot(distance,10*log10(signal_est_periodogram(:,1)));
% plot(distance,10*log10(signal_est_periodogram_correct(:,1)));
plot(distance, 10*log10(reference_spower));
plot(distance, ones(1,length(distance))*10*log10(tauc/pi));

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