%% written by Taehwan Kim (taehwan@berkeley.edu)
header;
c = 3e8;
% alpha = 10e6/(352/3e8);

% datapathname = '~/isg/MOABB/lidar_experiments/experiments/data/laser/linewidth/cleo_dec9/result_320m_100u_191mA_2.mat';
% datapathname = '~/isg/MOABB/lidar_experiments/experiments/data/laser/linewidth/ming_dec5/result_ming_300m_50u_40mA.mat';
% datapathname = strcat('~/isg/MOABB',dataname);
% addpath(datapathname);

%% open simlink model
open_system('./model/symodel_simple_lp');

%% set constants
q = 1.6e-19;
Res = 1;
c = 3e8;
timestep = 1e-9;

chirpbw = 15e9; %%%%%%%
realstop = 100*1e-6;
alpha = chirpbw/realstop;

distance = 200;
tau = 2*distance/c;
buffer = 5000;
                                                                                
stoptime = realstop + tau; %%%%%%%

Fs = 1/timestep;
omega0 = 0;
phi0 = 0;
prx=1/sqrt(2);

N = round(realstop/timestep);

% fnoise_list = 10.^(6:0.2:7); %%%%%%%
% fnoise_list = [5e5 1e6 5e6];
fnoise_list = 1e6;
% Prx = 1e-3;
% Prx = [10^-9 10^-10 10^-11];
% Prx = 10.^(-11:0.2:-8);
% Prx = 10^-10.6;
Prx_list = 1e-3;
% Prx_list = 10.^(-11:0.5:-8);
Plo = 10e-3;
% Prx = 1e-9;
% [10^-10 10^-9 10^-8];
% snoise_list = q/Res./Prx;
snoise = 2*q*Res*Plo;


%% run simulation

iter = 50;

periods = (1./(tau*alpha))/timestep;
rightlength = floor(periods.*floor(N./periods));

x_windowed = zeros(length(fnoise_list), iter, N);
signal_est_periodogram = zeros(length(fnoise_list), iter);
dist_est_periodogram = zeros(length(Prx_list),length(fnoise_list), iter);
dist_est_periodogram_lse = zeros(length(Prx_list),length(fnoise_list), iter);
dist_est_rb = zeros(length(Prx_list),length(fnoise_list), iter);

std_est_periodogram = zeros(length(Prx_list),length(fnoise_list));
std_est_periodogram_lse = zeros(length(Prx_list),length(fnoise_list));
std_est_rb = zeros(length(Prx_list),length(fnoise_list));
max_est_periodogram = zeros(length(Prx_list),length(fnoise_list));

percollection = zeros(length(Prx_list),length(fnoise_list),iter,N/2+1);
percollection_lp = zeros(length(Prx_list),length(fnoise_list),iter,N/2+1);
solution_lse_collection = zeros(length(Prx_list),length(fnoise_list),iter,2);

tolerance = 20;

rng(510);
A = randi(10000,2,iter);

for k=1:length(Prx_list)
    for j=1:length(fnoise_list)
        i = 1;
        fail_flag = 0;
        
        Prx = Prx_list(k);
        fnoise = fnoise_list(j);
        tauc=1./(pi*fnoise);
        dc = (c/(4*pi*fnoise))*log(realstop*pi*fnoise);
        
        while (i <= iter && fail_flag == 0)
            
%             snoise = snoise_list;
            fseed = A(1,i);
            sseed = A(2,i);
            
            
            outlier_flag = 1;
            outlier_count = 0;
            
            while (outlier_flag == 1)
                                
                sim('symodel_simple_lp');
                measdata_lp = squeeze(vout15.Data);
                measdata = squeeze(vout2.Data);
%                 noise_out = squeeze(vout.Data);
%                 lp_out = squeeze(vout12.Data);
                measdata = measdata(end-N:end-N+rightlength-1);
                measdata_lp = measdata_lp(end-N:end-N+rightlength-1);
                
                x_windowed(j,i,:)= cat(2,measdata',zeros(1,N-length(measdata)));
                x_windowed_lp(j,i,:)= cat(2,measdata_lp',zeros(1,N-length(measdata_lp)));

%                     zeros(round(4*N/5)-rightlength,1)
                tempx_windowed = squeeze(x_windowed(j,i,:));
                tempx_windowed_lp = squeeze(x_windowed_lp(j,i,:));
%                 [per_noise,Fn] = periodogram(noise_out,[],'onesided',length(noise_out),Fs);
%                 [per_lp,Fn] = periodogram(lp_out,[],'onesided',length(lp_out),Fs);

                [per_lp,F] = periodogram(tempx_windowed_lp,window(@rectwin,length(tempx_windowed_lp)),length(tempx_windowed_lp),Fs);
                [per,F] = periodogram(tempx_windowed,window(@rectwin,length(tempx_windowed)),length(tempx_windowed),Fs);

                [Y, I] = max(per);
                signal_est_periodogram(j,i) = Y;
                dist_est_periodogram(k,j,i) = F(I)*c/alpha/2;        
                
                dist_est_rb(k,j,i) = frequency_rb(tempx_windowed,Fs)*c/alpha/2;
                    
                solution_lse = pnoise_fit(10*log10(per),F,fnoise,snoise,[I 0.5*Y*(2*pi*fnoise)^2]);
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
            percollection(k,j,i,:) = per;
            percollection_lp(k,j,i,:) = per_lp;
            solution_lse_collection(k,j,i,:) = solution_lse;
            i = i + 1;
        end
        
        fprintf('Just finished iteration #%d-%d, fail_flag: %d, ', k,j,fail_flag);
        if fail_flag == 1
            std_est_periodogram_lse(k,j) = 100;
        else
            std_est_periodogram_lse(k,j) = std(dist_est_periodogram_lse(k,j,:));
        end
        max_est_periodogram(k,j) = mean(signal_est_periodogram(j,:));
        
        std_est_periodogram(k,j) = std(dist_est_periodogram(k,j,:));
        std_est_rb(k,j) = std(dist_est_rb(k,j,:));
        
        fprintf('lse: %.2f cm\n',single(std_est_periodogram_lse(k,j)*100));
        fprintf('per: %.2f cm\n',single(std_est_periodogram(k,j)*100));
        fprintf('rb: %.2f cm\n',single(std_est_rb(k,j)*100));
    end
end

figure();
hold on;
plot(F, 10*log10(squeeze(mean(percollection(1,1,:,:),3))));
plot(F, 10*log10(squeeze(mean(percollection_lp(1,1,:,:),3))));

% close all;
% fig1handle=fig('units','inches','width',2*textwidth,'height',textwidth*0.5,'font','Times','fontsize',2*fs);
% colormap(fig1handle,cm_viridis);
% subhandle1=subplot(1,3,1);
% set(gcf,'DefaultAxesColorOrder',colororder);
% hold on
% plot(F/1e6,30+10*log10(squeeze(percollection(5,1,1,:))),'color',cb);
% plot(F/1e6,30+10*log10(squeeze(mean(percollection(5,1,:,:),3))),'color',co);
% % plot(F/1e6,30+pnoise_model(fnoise,snoise,squeeze(solution_lse_collection(:,:,10,:)),F),'color',cg,'linewidth',2);
% % plot(F,pnoise_model(fnoise,snoise,[solution_lse(1) 2*Res^2*Plo*Prx*2*pi*fnoise],F));
% set(gca,'linewidth',1.5*lw);
% % set(gca,'yminortick','off');
% set(gca,'ticklength',[0.02 0.02]);
% xlabel('Frequency (MHz)');
% ylabel('Power Density (dBm/Hz)');
% % legend('Single trial','Averaged (50 trials)');
% legendstring={'Single trial','Averaged (50 trials)'};
% legendflex(legendstring,'anchor',{'ne','ne'},'buffer',[0 0],'bufferunit','inches','box','off','xscale',1);
% legend boxoff;
% box off;
% xlim([0 500]);
% ylim([-220 -100]);
% 
% subhandle2=subplot(1,3,2);
% set(gcf,'DefaultAxesColorOrder',colororder);
% edges=(198:0.1:202);
% hold on
% histogram(dist_est_periodogram(1,1,:),edges,'facecolor',cb);
% histogram(dist_est_periodogram_lse(1,1,:),edges,'facecolor',co);
% set(gca,'linewidth',1.5*lw);
% set(gca,'yminortick','off');
% set(gca,'ticklength',[0.02 0.02]);
% xlabel('Distance (m)');
% ylabel('Count');
% legendstring={'L-LSE','MF'};
% % legendflex(legendstring,'anchor',{'ne','ne'},'buffer',[0 0],'bufferunit','inches','box','off','xscale',0.5);
% legend(legendstring);
% legend boxoff;
% box off;
% xlim([198 202]);
% ylim([0 30]);
% 
% subhandle3 = subplot(1,3,3);
% mincolor = 0; maxcolor = 15;
% Prx = {'10^{-11}' '' '10^{-10}' '' '10^{-9}' '' '10^{-8}'};
% fnoise_list = {'10^{0}' '' '' '' '' '10^{1}'};
% hm1 = heatmapcus(flipud(single(std_est_periodogram_lse*100)),fnoise_list,fliplr(Prx), [],...
%         'ShowAllTicks', true,'TickTexInterpreter','true',...
%         'MinColorValue', mincolor, 'MaxColorValue', maxcolor,'GridLines','-');
% set(gca,'ticklength',[0 0]);
% colorbar;
% xlabel('Linewidth (MHz)')
% ylabel('RX Power (W)')
% 
% h1=get(subhandle1,'Position');
% h2=get(subhandle2,'Position');
% h3=get(subhandle3,'Position');
% 
% plotx1 = h1(1)-0.05;
% ploty1 = h1(2)+0.05;
% 
% plotwidth1 = h1(3)*1.15;
% plotheight1 = h1(4)*0.95;
% % plotpitch = plotwidth + 0.025;
% 
% 
% plotx2 = h2(1);
% ploty2 = h2(2)+0.05;
% 
% plotwidth2 = h2(3)*1.15;
% plotheight2 = h2(4)*0.95;
% % plotpitch = plotwidth + 0.025;
% 
% 
% plotx3 = h3(1)+0.05;
% ploty3 = h3(2)+0.05;
% 
% plotwidth3 = h3(3)*1.05;
% plotheight3 = h3(4)*0.95;
% % plotpitch = plotwidth + 0.025;
% 
% set(subhandle1,'Position',[plotx1 ploty1 plotwidth1 plotheight1]);
% set(subhandle2,'Position',[plotx2 ploty2 plotwidth2 plotheight2]);
% set(subhandle3,'Position',[plotx3 ploty3 plotwidth3 plotheight3]);
% 
% pos = [3.2 0.5];
% str = 'Accuracy (cm)';
% t=text(subhandle3,pos(1),pos(2), str,'fontname','times',...
%     'fontsize',2*fs,'units','inches','rotation',90);
% 
% pos = [-0.75 -0.5];
% str = '(a)';
% t=text(subhandle1,pos(1),pos(2), str,'fontname','times',...
%     'fontsize',2*fs,'units','inches');
% 
% pos = [-0.75 -0.5];
% str = '(b)';
% t=text(subhandle2,pos(1),pos(2), str,'fontname','times',...
%     'fontsize',2*fs,'units','inches');
% % 
% % pos = [0 1];
% % str = '\sigma_{LLSE}=7cm\n\sigma_{MF}=7cm';
% % t=text(subhandle2,pos(1),pos(2), str,'fontname','times',...
% %     'fontsize',2*fs,'units','inches');
% 
% pos = [-0.75 -0.5];
% str = '(c)';
% t=text(subhandle3,pos(1),pos(2), str,'fontname','times',...
%     'fontsize',2*fs,'units','inches');
% 
% boxw = 0.1;
% boxh = 0.25;
% dim = [0.42 0.55 boxw boxh];
% str = {'$\sigma_{LLSE}$=3.56cm'};
% annotation(...
%     'textbox',dim, 'String',str,'FitBoxToText','off','fontname','times',...
%     'fontsize',1.5*fs,'horizontalalignment','center','Color',ck,...
%     'verticalalignment','middle','linestyle','none','Interpreter','latex');
% 
% dim = [0.42 0.5 boxw boxh];
% str = {'$\sigma_{MF}$=37.9cm'};
% annotation(...
%     'textbox',dim, 'String',str,'FitBoxToText','off','fontname','times',...
%     'fontsize',1.5*fs,'horizontalalignment','center','Color',ck,...
%     'verticalalignment','middle','linestyle','none','Interpreter','latex');
% 
% 
% 
% 
% save2pdf('./cleo2018/fig1.pdf',fig1handle,600);

%%
% sample = load(datapathname);
% 
% % filename = sprintf('./data/pc_10m_50khz_133mvpp_820mv/scope_26.csv');
% wmeas = transpose(sample.wmeas(:,:));
% tmeas = transpose(sample.tmeas(2:end,:));
% 
% wstd = sqrt(var(wmeas,[],2));
% nor = (1/wstd(1))*eye(length(wstd));
% % nor = diag(1./wstd);
% wmeas = nor*wmeas;
% leng = length(tmeas)-1;
% wmeas = wmeas(:,1:leng);
% tmeas = tmeas(:,1:leng);
% % wmeas = wmeas(:,end-leng:end);
% % tmeas = tmeas(:,end-leng:end);
% % scope = csvread(filename,5,0);
% % scope = scope';
% [lenx, leny] = size(wmeas);
% 
% [per,F]=periodogram(wmeas(1,:),window(@rectwin,length(wmeas(1,:))),length(wmeas(1,:)),1/(tmeas(2)-tmeas(1)));
% 
% % data = zeros(lenx,leny);
% per_data = zeros(lenx,length(per)-49-50);
% % per_data_lp = zeros(lenx,length(per)-49-300);
% % window_lp=zeros(1,size(per_data,2));
% % template = normpdf(linspace(1,3898,3898),0,3);
% % window_lp(1:1949)= template(1:1949);
% % window_lp(end-1948:end)=fliplr(template(1:1949));
% 
% solutions_per = zeros(lenx,1);
% % solutions_per_lp = zeros(lenx,1);
% solutions_rb = zeros(lenx,1);
% solutions_lse = zeros(lenx,1);
% solution_lse = zeros(lenx,2);
% tau = zeros(lenx,1);
% 
% roisize=500;
% for i=1:lenx
% % for i=92
% 
% %     filename = sprintf('scope_%d.csv', i+25);
% %     filename = strcat(datapathname, filename);
% %     scope = csvread(filename,5,0);
% %     scope = scope';
% %     data(i,:) = wmeas(i,:);
%     [per_temp,f_temp] = periodogram(wmeas(i,:),window(@rectwin,length(wmeas(i,:))),length(wmeas(i,:)),1/(tmeas(2)-tmeas(1)));
%     per_temp = per_temp(50:end-50);
%     f_temp = f_temp(50:end-50);
%     per_data(i,:) = per_temp;
%     
% %     per_data_lp(i,:) = abs(ifft(window_lp.*fft(per_data(i,:))));
%     
% %     f_temp = (f_temp(1:2000));
% %     per_temp = per_temp(1:2000);
%     
%     [Y, I] = max(per_data(i,:));
%     
%     solutions_per(i) = f_temp(I);
% %     solutions_rb(i) = frequency_rb(wmeas(i,:),1/(tmeas(2)-tmeas(1)));
%     f_roi = f_temp(I-roisize:I+roisize-1);
%     per_temp_roi = per_temp(I-roisize:I+roisize-1);
%     
%     solution_lse(i,:) = pnoise_fit_data(10*log10(per_temp_roi),f_roi,[roisize Y]);
% %     solution_lse(i,:) = pnoise_fit_data(10*log10(per_temp),f_temp,[I Y (mean(per_data(end-30:end)))]);
% 
%     solutions_lse(i) = f_roi(1) + solution_lse(i,1)*(f_roi(2)-f_roi(1));
% %     tau(i) = solution_lse(3);
%     
% %     [Y, I] = max(per_data_lp(i,:));
% %     solutions_per_lp(i) = f_temp(I);
% end
% 
% alpha = mean(solutions_per)/(340*10/8/3e8);
% 
% alpha = 1.279e7/(304*1.49/3e8);

% averaged_per = mean(per_data,1);
% pers = zeros(9,size(averaged_per,2));
% pers(9,:) = averaged_per;
%     
% figure(1);
% hold on;
% % plot(f_temp,10*log10(pnoise_model_data(solution_lse(end,:),f_temp)));
% plot(f_temp,10*log10(averaged_per));
% plot(f_roi,(pnoise_model_data(solution_lse(1,:),f_roi)));
% plot(f_roi,(pnoise_model_data(solution_lse(10,:),f_roi)));
% plot(f_roi,(pnoise_model_data(solution_lse(20,:),f_roi)));
% std(solutions_per)/std(solutions_lse)
% % xlim([0 2e7]);
% 
% % mean(solutions_per)*c/alpha/2
% % mean(solutions_lse)*c/alpha/2
% % 
% % d_lse = solutions_lse*c/alpha/2;
% std(solutions_lse)*c/alpha/2
% std(solutions_per)*c/alpha/2
% % % 
% edges=(200:1:240);
% figure(2);
% hold on;
% histogram(solutions_per*c/alpha/2,edges);
% histogram(solutions_lse*c/alpha/2,edges);
% histogram(solutions_per_lp*c/alpha/2,edges);

% plot(F,10*log10(averaged_per))
% 
% legend('30m','90m','150m','210m','270m');
% % F = F(1:2000);
% F = F';
% averaged_per = averaged_per(1:length(F));
% 
% solution_lse = pnoise_fit_data(averaged_per,F,[1000 0 0.7*8.626460257595531e-09]);
% solution = pnoise_model_data(solution_lse,F);

% figure(1);

