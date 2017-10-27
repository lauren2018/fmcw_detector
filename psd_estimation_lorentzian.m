%% written by Taehwan Kim (taehwan@berkeley.edu)
clear; close all; clc;
format long;
set(0,'DefaultFigureWindowStyle','normal')
addpath './frequency'
addpath './Colormaps'
addpath './distinguishable_colors/'
addpath './hex_and_rgb_v2'
addpath './fig'
addpath './legendflex'
addpath './setgetpos_V1.2'

header;

m=100;
cm_magma=magma(m);
cm_inferno=inferno(m);
cm_plasma=plasma(m);
cm_viridis=viridis(m);

%% open simulink model
open_system('./model/symodel_simple');



%% set constants
q = 1.6e-19;
Res = 1;
c = 3e8;
timestep = 0.5e-9;
stoptime = 10e-6; %%%%%%%
Fs = 1/timestep;
chirpbw = 10e9; %%%%%%%
alpha = chirpbw/stoptime;
omega0 = 0;
phi0 = 0;
prx=1/sqrt(2);

N = round(stoptime/timestep);
freq = 0:Fs/N:Fs/2;
fnoise = 25e6; %%%%%%%
tauc=1/(pi*fnoise);
dc = (c/(4*pi*fnoise))*log(stoptime*pi*fnoise);
% Prx = 10.^(-11:0.2:-9);
% Prx = 0.316227766016840e-9;
% Prx = [10^-3 10^-9];
Prx = 10^-3;
snoise_list = q/Res./Prx;
% snoise_list =0;
% snoise = 0;
amplitude=1;

buffer = 4096;

%% run simulation
% distance = 20;
% distance = 20;
iter = 100;
padm = 0;

% distance=100;
% distance = 10.1:10:100.1;
% distance = [11 22 33 44 55 66 77 88 99 110];
distance = 120;

dist_res = 0.5*c/chirpbw;
% distance = [3 4];
tau_list = 2*distance/c;
reference_spower = stoptime * exp(-2*tau_list/tauc) / pi;

periods = (1./(tau_list*alpha))/timestep;
rightlength = floor(periods.*floor(0.9*N./periods));

x = zeros(length(distance), iter, (1+2*padm)*N);
x_windowed = zeros(length(distance), iter, (1+2*padm)*N);
signal_est_periodogram = zeros(length(distance), iter);

dist_est_periodogram = zeros(length(Prx),length(distance), iter);
dist_est_periodogram_lse = zeros(length(Prx),length(distance), iter);
dist_est_periodogram_mle = zeros(length(Prx),length(distance), iter);
dist_est_music = zeros(length(Prx),length(distance), iter);
dist_est_periodogram_correct = zeros(length(Prx),length(distance), iter);
dist_est_rb = zeros(length(Prx),length(distance), iter);

std_est_periodogram = zeros(length(Prx),length(distance));
std_est_periodogram_lse = zeros(length(Prx),length(distance));
std_est_periodogram_mle = zeros(length(Prx),length(distance));
std_est_music = zeros(length(Prx),length(distance));
std_est_rb = zeros(length(Prx),length(distance));

max_est_periodogram = zeros(length(Prx),length(distance));

percollection = zeros(length(Prx),length(distance),iter,N/2+1);

tolerance = 5;

for k=1:length(Prx)
    for j=1:length(distance)
        for i=1:iter
            snoise = snoise_list(k);
            tau = tau_list(j);
            fseed = 3*i+1;
            sseed = 5*i+1;
            outlier_flag = 1;
            
            while (outlier_flag == 1)
                sim('symodel_simple');
                measdata = squeeze(vout15.Data);
                x_windowed(j,i,:)=padarray(cat(1,zeros(round(N/10),1),measdata(round(N/10)+1:round(N/10)+rightlength(j)),zeros(round(9*N/10)-rightlength(j),1)),padm*N,0,'both');

                tempx_windowed = squeeze(x_windowed(j,i,:));

                [per,F] = periodogram(tempx_windowed,window(@rectwin,length(tempx_windowed)),length(tempx_windowed),Fs);
                
                

                [Y, I] = max(per);
                signal_est_periodogram(j,i) = Y;
                dist_est_periodogram(k,j,i) = F(I)*c/alpha/2;        

                solution_lse = pnoise_fit(per,F,[N/4 tauc snoise]);
                dist_est_periodogram_lse(k,j,i) = c*solution_lse(1)*(F(2)-F(1))/alpha/2;

                solution_mle = mlelorentz(per,F,[N/4 tauc snoise]);
                dist_est_periodogram_mle(k,j,i) = c*solution_mle(1)*(F(2)-F(1))/alpha/2;

                dist_est_rb(k,j,i) = frequency_rb(tempx_windowed,Fs)*c/alpha/2;

                [S, F] = pmusic(tempx_windowed,6,freq,Fs);
                [Y, I] = max(S);
                dist_est_music(k,j,i) = F(I)*c/alpha/2;
                
                outlier_flag_periodogram = outlier_test(distance, dist_est_periodogram(k,j,i), tolerance);
                outlier_flag_music = outlier_test(distance, dist_est_music(k,j,i), tolerance);
                outlier_flag_lse = outlier_test(distance, dist_est_periodogram_lse(k,j,i), tolerance);
                
                if (outlier_flag_lse + outlier_flag_music + outlier_flag_periodogram == 0)
                    outlier_flag = 0;
                    fprintf('No outliers! %d-%d-%d ',outlier_flag_lse,outlier_flag_music,outlier_flag_periodogram);
                else
                    fprintf('Outlier exists %d-%d-%d\n',outlier_flag_lse,outlier_flag_music,outlier_flag_periodogram);
                    fseed = fseed + 1;
                    sseed = sseed + 1;
                end
            end
            fprintf('Just finished iteration #%d-%d-%d\n', k,j,i);
            
            percollection(k,j,i,:) = per;
            
        end
        std_est_periodogram(k,j) = std(dist_est_periodogram(k,j,:));
        std_est_periodogram_lse(k,j) = std(dist_est_periodogram_lse(k,j,:));
        std_est_periodogram_mle(k,j) = std(dist_est_periodogram_mle(k,j,:));
        std_est_music(k,j) = std(dist_est_music(k,j,:));
        std_est_rb(k,j) = std(dist_est_rb(k,j,:));
        max_est_periodogram(k,j) = mean(signal_est_periodogram(j,:));
    end
end
% save('data.mat','x');
    
%% plot the periodogram-estmated psd for 

temp = squeeze(percollection);
temp = mean(temp,1);
figure
hold on;
plot(F,10*log10(temp));
solution_lse = pnoise_fit(temp',F,[N/4 tauc snoise]);
solution_lse(2) = solution_lse(2)*4;
plot(F,-3+10*log10(pnoise_model(solution_lse,F)));

