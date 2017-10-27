%% written by Taehwan Kim (taehwan@berkeley.edu)
header;

dataname = 'pc_10m_50khz_133mvpp_820mv/';
datapathname = strcat('./data/',dataname);
% addpath(datapathname);

c = 3e8;
alpha = 100e9/10e-6;

filename = sprintf('./data/pc_10m_50khz_133mvpp_820mv/scope_26.csv');
scope = csvread(filename,5,0);
scope = scope';
[lenx, leny] = size(scope);
time = scope(1,:);

[per,F]=periodogram(scope(2,:),window(@rectwin,length(scope(2,:))),length(scope(2,:)),1/(time(2)-time(1)));

data = zeros(30,leny);
per_data = zeros(30,length(per));

solutions_per = zeros(30,1);
solutions_lse = zeros(30,1);

for i=1:30

    filename = sprintf('scope_%d.csv', i+25);
    filename = strcat(datapathname, filename);
    scope = csvread(filename,5,0);
    scope = scope';
    data(i,:) = scope(2,:);
    [per_temp,f_temp] = periodogram(data(i,:),window(@rectwin,length(data(i,:))),length(data(i,:)),1/(time(2)-time(1)));
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


