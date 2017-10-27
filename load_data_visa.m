%% written by Taehwan Kim (taehwan@berkeley.edu)
header;
c = 3e8;
alpha = 100e9/10e-6;

datapathname = '~/isg/MOABB/lidar_experiments/experiments/data/laser/linewidth/laser_beat_no_current_dbr/result_2017-10-24-17-41.mat';
% datapathname = strcat('~/isg/MOABB',dataname);
% addpath(datapathname);
sample = load(datapathname);

% filename = sprintf('./data/pc_10m_50khz_133mvpp_820mv/scope_26.csv');
wmeas = sample.wmeas';
tmeas = sample.tmeas';

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
    
    f_temp = (f_temp(1:2000));
    per_temp = per_temp(1:2000);
    
    [Y, I] = max(per_temp);
    solutions_per(i) = f_temp(I)*c/alpha;
    
    solution_lse = pnoise_fit_data(per_temp,f_temp,[1000 0]);
    solutions_lse(i) = c*solution_lse(1)*(f_temp(2)-f_temp(1))/alpha;
    
end

averaged_per = mean(per_data,1);

F = F(1:2000);
F = F';
averaged_per = averaged_per(1:2000);

solution_lse = pnoise_fit_data(averaged_per,F,[1000 0]);
solution = pnoise_model_data(solution_lse,F);

figure();
hold on;
plot(F,10*log10(averaged_per));
plot(F,10*log10(solution));

fprintf('PER: %f,%f\n', mean(solutions_per), std(solutions_per));
fprintf('LLSE: %f, %f\n', mean(solutions_lse),std(solutions_lse));


