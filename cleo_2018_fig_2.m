%% written by Taehwan Kim (taehwan@berkeley.edu)
header;
c = 3e8;
% alpha = 10e6/(352/3e8);

datapathname = '~/isg/MOABB/lidar_experiments/experiments/data/laser/linewidth/cleo_dec11/result_12.mat';
% datapathname = '~/isg/MOABB/lidar_experiments/experiments/data/laser/linewidth/ming_dec5/result_ming_300m_50u_40mA.mat';
% datapathname = strcat('~/isg/MOABB',dataname);
% addpath(datapathname);
sample = load(datapathname);

% filename = sprintf('./data/pc_10m_50khz_133mvpp_820mv/scope_26.csv');
wmeas = transpose(sample.wmeas(:,:));
tmeas = transpose(sample.tmeas(2:end,:));

wstd = sqrt(var(wmeas,[],2));
nor = (1/wstd(1))*eye(length(wstd));
% nor = diag(1./wstd);
wmeas = nor*wmeas;
leng = length(tmeas)-1;
wmeas = wmeas(:,1:leng);
tmeas = tmeas(:,1:leng);
% wmeas = wmeas(:,end-leng:end);
% tmeas = tmeas(:,end-leng:end);
% scope = csvread(filename,5,0);
% scope = scope';
[lenx, leny] = size(wmeas);

[per,F]=periodogram(wmeas(1,:),window(@rectwin,length(wmeas(1,:))),length(wmeas(1,:)),1/(tmeas(2)-tmeas(1)));

% data = zeros(lenx,leny);
per_data = zeros(length(per)-199-50,lenx);
roisize=1500;
f_roi = zeros(2*roisize,lenx);
per_data_roi = zeros(2*roisize,lenx);
I=zeros(1,lenx);
% per_data_lp = zeros(lenx,length(per)-49-300);
% window_lp=zeros(1,size(per_data,2));
% template = normpdf(linspace(1,3898,3898),0,3);
% window_lp(1:1949)= template(1:1949);
% window_lp(end-1948:end)=fliplr(template(1:1949));

solutions_per = zeros(lenx,1);
% solutions_per_lp = zeros(lenx,1);
solutions_rb = zeros(lenx,1);
solutions_lse = zeros(lenx,1);
solution_lse = zeros(lenx,2);
tau = zeros(lenx,1);

for i=1:lenx
% for i=92

%     filename = sprintf('scope_%d.csv', i+25);
%     filename = strcat(datapathname, filename);
%     scope = csvread(filename,5,0);
%     scope = scope';
%     data(i,:) = wmeas(i,:);
    [per_temp,f_temp] = periodogram(wmeas(i,:),window(@rectwin,length(wmeas(i,:))),length(wmeas(i,:)),1/(tmeas(2)-tmeas(1)));
    per_temp = per_temp(200:end-50);
    f_temp = f_temp(200:end-50);
    per_data(:,i) = per_temp;
    
%     per_data_lp(i,:) = abs(ifft(window_lp.*fft(per_data(i,:))));
    
%     f_temp = (f_temp(1:2000));
%     per_temp = per_temp(1:2000);
    
    [Y, I(i)] = max(per_data(:,i));
    
    solutions_per(i) = f_temp(I(i));
% %     solutions_rb(i) = frequency_rb(wmeas(i,:),1/(tmeas(2)-tmeas(1)));


    f_roi(:,i) = f_temp(I(i)-roisize:I(i)+roisize-1);
    per_temp_roi(:,i) = per_temp(I(i)-roisize:I(i)+roisize-1);
    

    solution_lse(i,:) = pnoise_fit_data((per_temp_roi(:,i)),f_roi(:,i),[roisize 0.1*Y]);


% %     solution_lse(i,:) = pnoise_fit_data(10*log10(per_temp),f_temp,[I Y (mean(per_data(end-30:end)))]);
% f_roi(1) + 


    solutions_lse(i) = f_roi(1,i)+solution_lse(i,1)*(f_roi(2,i)-f_roi(1,i));
%     solutions_lse(i) = solution_lse(i,1)*(f_temp(2)-f_temp(1));


% %     tau(i) = solution_lse(3);
%     
% %     [Y, I] = max(per_data_lp(i,:));
% %     solutions_per_lp(i) = f_temp(I);
end

alpha = mean(solutions_per)/(220/3e8);

averaged_per = mean(per_data,2);



close all;
fig1handle=fig('units','inches','width',2*textwidth,'height',textwidth*0.4,'font','Times','fontsize',2*fs);
colormap(fig1handle,cm_viridis);

subhandle1=subplot(1,3,1);
fgimage=imread('./setup.png');
imshow(fgimage);
% hfg.AlphaData = alpha;
axis image;
% set(gca,'backgroundcolor',cb)

subhandle2=subplot(1,3,2);
set(gcf,'DefaultAxesColorOrder',colororder);
hold on
% f_temp,10*log10(averaged_per)
plot(f_temp/1e6,30+10*log10(per_data(:,i)),'color',cb);
% plot(f_temp/1e6,30+10*log10(averaged_per),'color',cy);
% plot(f_roi(:,:)/1e6,30+10*log10(,'color',co);
plot(f_roi(:,i)/1e6,30+10*log10(pnoise_model_data(solution_lse(i,:),f_roi(:,i))),'color',co,'linewidth',2);

% plot(F/1e6,30+pnoise_model(fnoise,snoise,squeeze(solution_lse_collection(:,:,10,:)),F),'color',cg,'linewidth',2);
% plot(F,pnoise_model(fnoise,snoise,[solution_lse(1) 2*Res^2*Plo*Prx*2*pi*fnoise],F));
set(gca,'linewidth',1.5*lw);
% set(gca,'yminortick','off');
set(gca,'ticklength',[0.02 0.02]);
xlabel('Frequency (MHz)');
ylabel('PSD (dBm/Hz)');
% legend('Single trial','Averaged (50 trials)');
legendstring={'Single trial','Fitted'};
legendflex(legendstring,'anchor',{'ne','ne'},'buffer',[0 0],'bufferunit','inches','box','off','xscale',1,'fontsize',1.5*fs);
legend boxoff;
box off;
xlim([50 300]);
ylim([-80 -20]);

subhandle3=subplot(1,3,3);
set(gcf,'DefaultAxesColorOrder',colororder);
% edges=(198:0.1:202);
edges=(90:1:130);
hold on
% histogram(solutions_per*c/alpha/2,edges);
% histogram(solutions_lse*c/alpha/2,edges);

histogram(solutions_lse*c/alpha/2,edges,'facecolor',cb);
histogram(solutions_per*c/alpha/2,edges,'facecolor',co);
set(gca,'linewidth',1.5*lw);
set(gca,'yminortick','off');
set(gca,'ticklength',[0.02 0.02]);
xlabel('Distance (m)');
ylabel('Count');
legendstring={'L-LSE','MF'};
% legendflex(legendstring,'anchor',{'ne','ne'},'buffer',[0 0],'bufferunit','inches','box','off','xscale',0.5);
legend(legendstring);
legend boxoff;
box off;
xlim([edges(1) edges(end)]);
ylim([0 40]);

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

h1=get(subhandle1,'Position');
h2=get(subhandle2,'Position');
h3=get(subhandle3,'Position');

plotx1 = h1(1)-0.1;
ploty1 = h1(2)-0.25;

plotwidth1 = h1(3)*1.5;
plotheight1 = h1(4)*1.5;
% plotpitch = plotwidth + 0.025;


plotx2 = h2(1);
ploty2 = h2(2)+0.05;

plotwidth2 = h2(3)*1.15;
plotheight2 = h2(4)*0.95;
% plotpitch = plotwidth + 0.025;


plotx3 = h3(1)+0.05;
ploty3 = h3(2)+0.05;

plotwidth3 = h3(3)*1.15;
plotheight3 = h3(4)*0.95;
% plotpitch = plotwidth + 0.025;

set(subhandle1,'Position',[plotx1 ploty1 plotwidth1 plotheight1]);
set(subhandle2,'Position',[plotx2 ploty2 plotwidth2 plotheight2]);
set(subhandle3,'Position',[plotx3 ploty3 plotwidth3 plotheight3]);

% pos = [3.2 0.5];
% str = 'Accuracy (cm)';
% t=text(subhandle3,pos(1),pos(2), str,'fontname','times',...
%     'fontsize',2*fs,'units','inches','rotation',90);

% pos = [-0.75 -0.5];
% str = '(a)';
% t=text(subhandle1,pos(1),pos(2), str,'fontname','times',...
%     'fontsize',2*fs,'units','inches');
% 
% pos = [-0.75 -0.5];
% str = '(b)';
% t=text(subhandle2,pos(1),pos(2), str,'fontname','times',...
%     'fontsize',2*fs,'units','inches');
% 
% pos = [0 1];
% str = '\sigma_{LLSE}=7cm\n\sigma_{MF}=7cm';
% t=text(subhandle2,pos(1),pos(2), str,'fontname','times',...
%     'fontsize',2*fs,'units','inches');

% pos = [-0.75 -0.5];
% str = '(c)';
% t=text(subhandle3,pos(1),pos(2), str,'fontname','times',...
%     'fontsize',2*fs,'units','inches');


boxw = 0.1;
boxh = 0.1;
dim = [0.02 0.01 boxw boxh];
str = {'(a)'};
annotation(...
    'textbox',dim, 'String',str,'FitBoxToText','off','fontname','times',...
    'fontsize',2*fs,'horizontalalignment','left','Color',ck,...
    'verticalalignment','middle','linestyle','none');

boxw = 0.1;
boxh = 0.1;
dim = [0.35 0.01 boxw boxh];
str = {'(b)'};
annotation(...
    'textbox',dim, 'String',str,'FitBoxToText','off','fontname','times',...
    'fontsize',2*fs,'horizontalalignment','left','Color',ck,...
    'verticalalignment','middle','linestyle','none');

boxw = 0.1;
boxh = 0.1;
dim = [0.68 0.01 boxw boxh];
str = {'(c)'};
annotation(...
    'textbox',dim, 'String',str,'FitBoxToText','off','fontname','times',...
    'fontsize',2*fs,'horizontalalignment','left','Color',ck,...
    'verticalalignment','middle','linestyle','none');

boxw = 0.1;
boxh = 0.25;
dim = [0.75 0.6 boxw boxh];
str = {'$\sigma_{LLSE}$=43.9cm'};
annotation(...
    'textbox',dim, 'String',str,'FitBoxToText','off','fontname','times',...
    'fontsize',1.5*fs,'horizontalalignment','center','Color',ck,...
    'verticalalignment','middle','linestyle','none','Interpreter','latex');

dim = [0.75 0.5 boxw boxh];
str = {'$\sigma_{MF}$=2.9m'};
annotation(...
    'textbox',dim, 'String',str,'FitBoxToText','off','fontname','times',...
    'fontsize',1.5*fs,'horizontalalignment','center','Color',ck,...
    'verticalalignment','middle','linestyle','none','Interpreter','latex');




% save2pdf('./cleo2018/fig1.pdf',fig1handle,600);

% pers = zeros(9,size(averaged_per,2));
% pers(9,:) = averaged_per;
    
% solutions_lse = solutions_lse(abs(3e8*solutions_lse/alpha-220)<10);


figure(1);
hold on;
% plot(f_temp,10*log10(pnoise_model_data(solution_lse(end,:),f_temp)));
plot(f_temp,10*log10(per_data(:,i)));
plot(f_roi(:,i),10*log10(per_temp_roi(:,i)));
plot(f_temp,10*log10(averaged_per));
% plot(f_roi(:,end-2),10*log10(per_temp_roi(:,end-2)));
% plot(f_temp,10*log10(pnoise_model_data(solution_lse(1,:),f_temp)));
plot(f_roi(:,i),10*log10(pnoise_model_data(solution_lse(i,:),f_roi(:,i))));
% plot(f_roi,(pnoise_model_data(solution_lse(20,:),f_roi)));
std(solutions_per)/std(solutions_lse)
% xlim([0 2e7]);
% 
% % mean(solutions_per)*c/alpha/2
% % mean(solutions_lse)*c/alpha/2
% % 
% % d_lse = solutions_lse*c/alpha/2;
% std(solutions_lse)*c/alpha/2
% std(solutions_per)*c/alpha/2
% % % 
% edges=(80:1:150);
% figure(2);
% hold on;
% histogram(solutions_per*c/alpha/2,edges);
% histogram(solutions_lse*c/alpha/2,edges);
% % histogram(solutions_per_lp*c/alpha/2,edges);

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

