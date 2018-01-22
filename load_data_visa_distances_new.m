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

% mean(solutions_per)*c/alpha/2
% mean(solutions_lse)*c/alpha/2
% 
% d_lse = solutions_lse*c/alpha/2;
std(solutions_lse)*c/alpha/2
std(solutions_per)*c/alpha/2
% % 
edges=(80:1:150);
figure(2);
hold on;
histogram(solutions_per*c/alpha/2,edges);
histogram(solutions_lse*c/alpha/2,edges);
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

