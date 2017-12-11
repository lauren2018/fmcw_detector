%% written by Taehwan Kim (taehwan@berkeley.edu)
header;
c = 3e8;
alpha = 100e9/10e-6;

datapathname = '~/isg/MOABB/lidar_experiments/experiments/data/laser/linewidth/lidar_1102_8mA_133vpp_6250_10us_22p5celcius/result_4p5m_resampled.mat';
% datapathname = strcat('~/isg/MOABB',dataname);
% addpath(datapathname);
sample = load(datapathname);

% filename = sprintf('./data/pc_10m_50khz_133mvpp_820mv/scope_26.csv');
wmeas = sample.wmeas;
tmeas = sample.tmeas;

% scope = csvread(filename,5,0);
% scope = scope';
[lenx, leny] = size(wmeas);

[per,F]=periodogram(wmeas(1,:),window(@rectwin,length(wmeas(1,:))),length(wmeas(1,:)),1/(tmeas(2)-tmeas(1)));

% data = zeros(lenx,leny);
per_data = zeros(lenx,length(per));

solutions_per = zeros(lenx,1);
solutions_lse = zeros(lenx,1);

for i=1:lenx

%     filename = sprintf('scope_%d.csv', i+25);
%     filename = strcat(datapathname, filename);
%     scope = csvread(filename,5,0);
%     scope = scope';
%     data(i,:) = wmeas(i,:);
    [per_temp,f_temp] = periodogram(wmeas(i,:),window(@rectwin,length(wmeas(i,:))),length(wmeas(i,:)),1/(tmeas(2)-tmeas(1)));
    per_data(i,:) = per_temp;
    
%     f_temp = (f_temp(1:2000));
%     per_temp = per_temp(1:2000);
    
    [Y, I] = max(per_temp);
    solutions_per(i) = f_temp(I)*c/alpha/2;
    
    solution_lse = pnoise_fit_data(per_temp,f_temp,[1000 0 0.7*8.626460257595531e-09]);
    solutions_lse(i) = c*solution_lse(1)*(f_temp(2)-f_temp(1))/alpha/2;
    
end

averaged_per = mean(per_data,1);

% F = F(1:2000);
F = F';
averaged_per = averaged_per(1:length(F));

solution_lse = pnoise_fit_data(averaged_per,F,[1000 0 0.7*8.626460257595531e-09]);
solution = pnoise_model_data(solution_lse,F);

% figure(1);


fig1handle=fig('units','inches','width',textwidth/2,'height',textwidth/3,'font','Times','fontsize',fs);
set(gcf,'DefaultAxesColorOrder',colororder);
% figure(2);
hold on;
% plot(F/1e9,10*log10(per_temp));
plot(F/1e9,10*log10(averaged_per));
plot(F/1e9,10*log10(solution));
% 10*log10

legend('Single trial','Averaged','Fitted');
legend boxoff;
legend('Averaged','Fitted');
legend boxoff;
box off;
xlabel('Frequency (GHz)');
ylabel('Power density (W/Hz)');

% lbl = {'LLSE','R&B','MUSIC'};
% legendflex(lbl,'anchor',{'nw','nw'},'buffer',[0.2 0],'box','off','xscale',0.5);

set(gca,'linewidth',lw);
set(gca,'yminortick','off');
set(gca,'ticklength',[0.02 0.02]);
axis([0 0.3 -100 -70]);
% xticks([0 50 100]);
% yticks([0 0.5 1]);

% ylabel('count');


fprintf('PER: %f,%f\n', mean(solutions_per), std(solutions_per));
fprintf('LLSE: %f, %f\n', mean(solutions_lse),std(solutions_lse));
% save2pdf('psd_4p5m_avg.pdf',fig1handle,600);

fig2handle=fig('units','inches','width',textwidth/2,'height',textwidth/3,'font','Times','fontsize',fs);
set(gcf,'DefaultAxesColorOrder',colororder);

edges = (0.1:0.04:5);
hold on;
histogram(solutions_lse,edges);
histogram(solutions_per,edges);
axis([0 5 0 100]);
legend('LLSE','R&B');
legend boxoff;

xlabel('Distance (m)');
ylabel('Count');
% save2pdf('distance_4p5m.pdf',fig2handle,600);
