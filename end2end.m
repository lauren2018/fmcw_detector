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
timestep = 0.2e-9;

chirpbw = 10e9; %%%%%%%
realstop = 10*1e-6;
alpha = chirpbw/realstop;

distance = 75;
tau = 2*distance/c;
buffer = 14336;
                                                                                
stoptime = realstop + tau; %%%%%%%

Fs = 1/timestep;
omega0 = 0;
phi0 = 0;
prx=1/sqrt(2);

N = round(realstop/timestep);

% fnoise_list = 10.^(6:0.1:7); %%%%%%%
% fnoise_list = [5e5 1e6 5e6];
fnoise_list = 1e6;
% Prx = 1e-3;
% Prx = [10^-9 10^-10 10^-11];
% Prx = 10.^(-11:0.2:-8);
Prx = 10^-10.6;
% Prx = 1e-9;
% Prx = 10.^(-11.5:0.5:);
% Prx = 10^-3;
% [10^-10 10^-9 10^-8];
snoise_list = q/Res./Prx;


%% run simulation

iter = 50;

periods = (1./(tau*alpha))/timestep;
rightlength = floor(periods.*floor(N./periods));

x_windowed = zeros(length(fnoise_list), iter, N);
signal_est_periodogram = zeros(length(fnoise_list), iter);
dist_est_periodogram = zeros(length(Prx),length(fnoise_list), iter);
dist_est_periodogram_lse = zeros(length(Prx),length(fnoise_list), iter);
dist_est_rb = zeros(length(Prx),length(fnoise_list), iter);

std_est_periodogram = zeros(length(Prx),length(fnoise_list));
std_est_periodogram_lse = zeros(length(Prx),length(fnoise_list));
std_est_rb = zeros(length(Prx),length(fnoise_list));
max_est_periodogram = zeros(length(Prx),length(fnoise_list));

percollection = zeros(length(Prx),length(fnoise_list),N/2+1);

tolerance = 1;

rng(510);
A = randi(10000,2,iter);

for k=1:length(Prx)
    for j=1:length(fnoise_list)
        i = 1;
        fail_flag = 0;
        
        fnoise = fnoise_list(j);
        tauc=1./(pi*fnoise);
        dc = (c/(4*pi*fnoise))*log(realstop*pi*fnoise);
        
        while (i <= iter && fail_flag == 0)
            
            snoise = snoise_list(k);
            fseed = A(1,i);
            sseed = A(2,i);
            
            
            outlier_flag = 1;
            outlier_count = 0;
            
            while (outlier_flag == 1)
                                
                sim('symodel_simple');
                measdata = squeeze(vout15.Data);
                measdata = measdata(end-N:end-N+rightlength-1);
                
                x_windowed(j,i,:)= cat(2,measdata',zeros(1,N-length(measdata)));

%                     zeros(round(4*N/5)-rightlength,1)
                tempx_windowed = squeeze(x_windowed(j,i,:));

                [per,F] = periodogram(tempx_windowed,window(@rectwin,length(tempx_windowed)),length(tempx_windowed),Fs);

                [Y, I] = max(per);
                signal_est_periodogram(j,i) = Y;
                dist_est_periodogram(k,j,i) = F(I)*c/alpha/2;        
                
                dist_est_rb(k,j,i) = frequency_rb(tempx_windowed,Fs)*c/alpha/2;
                    
                solution_lse = pnoise_fit(per,F,[N/4 tauc snoise]);
                dist_est_periodogram_lse(k,j,i) = c*solution_lse(1)*(F(2)-F(1))/alpha/2;
                
                outlier_flag_lse = outlier_test(distance, dist_est_periodogram_lse(k,j,i), tolerance);
                
                if (outlier_flag_lse == 0)
                    outlier_flag = 0;
%                     fprintf('No outliers! ');
                elseif (outlier_count == 50)
                    fprintf('Failed! %d-%d\n', k,j);
                    outlier_flag = 0;
                    fail_flag = 1;
                else
                    fprintf('Outlier exists %d, %.2f\n',outlier_count,100*abs(distance-dist_est_periodogram_lse(k,j,i)));
                    fseed = fseed + 1;
                    sseed = sseed + 1;
                    outlier_count = outlier_count + 1;
                end
                
            end
            i = i + 1;
        end
        
        fprintf('Just finished iteration #%d-%d, fail_flag: %d, ', k,j,fail_flag);
        percollection(k,j,:) = per;
        if fail_flag == 1
            std_est_periodogram_lse(k,j) = 100;
        else
            std_est_periodogram_lse(k,j) = std(dist_est_periodogram_lse(k,j,:));
        end
        max_est_periodogram(k,j) = mean(signal_est_periodogram(j,:));
        
        std_est_periodogram(k,j) = std(dist_est_periodogram(k,j,:));
        std_est_rb(k,j) = std(dist_est_rb(k,j,:));
        
        fprintf('lse: %.2f cm\n',single(std_est_periodogram_lse(k,j)*100));
    end
end

hm=heatmap((fnoise_list), fliplr(Prx),flipud(single(std_est_periodogram_lse*100)));
colormap(hm,cm_magma);
caxis([0 20]);

