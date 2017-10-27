%% written by Taehwan Kim (taehwan@berkeley.edu)
header;
c = 3e8;
alpha = 100e9/10e-6;

datapathname = '~/isg/MOABB/lidar_experiments/experiments/data/laser/linewidth/result_2017-10-26-15-15.mat';
% datapathname = strcat('~/isg/MOABB',dataname);
% addpath(datapathname);
sample = load(datapathname);

% filename = sprintf('./data/pc_10m_50khz_133mvpp_820mv/scope_26.csv');
wmeas = sample.wmeas';
tmeas = sample.tmeas';

% scope = csvread(filename,5,0);
% scope = scope';
[lenx, leny] = size(wmeas);

tt = zeros(lenx,64);
for j=1:lenx
    wmt = wmeas(j,:);
    crossingtime = zeros(1,100);
    temp = zeros(1,leny);
    count = 1;

    for i=1:leny-1
        if wmt(i) * wmt(i+1) < 0
    %         if count ~= 0
    %         if tmeas(i) - crossingtime(count) > 0 
            count = count + 1;
            crossingtime(count) = tmeas(i);
            temp(i) = 1;
    %         end
    %         end
        end
    end

    crossingtime = crossingtime(1:count);
    periodsraw = crossingtime(2:end)-crossingtime(1:end-1);

    tempt = [];
    for i=1:length(periodsraw)
        if periodsraw(i) > 0.5e-7
            tempt = [tempt periodsraw(i)];
        end
    end

    tt(j,:) = tempt(2:65);
    j
%     lent = length(tt);
end
% realcrossing = zeros(1,cunt);
% count = 1;

% for i=1:length(crossingtime)-1
%     if crossingtime(i+1)-crossingtime(i) > 100*(tmeas(2)-tmeas(1))
%         realcrossing(count) = crossingtime(i+1);
%         count = count + 1;
%     end
% end
% 
% realcrossing = realcrossing(1:count-1);