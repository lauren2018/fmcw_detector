%% written by Taehwan Kim (taehwan@berkeley.edu)
header;
c = 3e8;
alpha = 100e9/10e-6;

datapathname = '~/isg/MOABB/lidar_experiments/experiments/data/laser/linewidth/result_2017-11-02-12-54.mat';
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

    crossingtime = crossingtime(2:count);
    periodsraw = crossingtime(2:end)-crossingtime(1:end-1);

    tempt = [];
    for i=1:length(periodsraw)
        if periodsraw(i) > 0.5e-7
            tempt = [tempt periodsraw(i)];
        end
    end

    tt(j,:) = tempt(2:65);
%     lent = length(tt);
end

tt_avg = mean(tt,1);

crosstimes = zeros(2,length(tt_avg));

for i=1:length(crosstimes)
    crosstimes(1,i) = sum(tt_avg(1:i));
    crosstimes(2,i) = i-1;
end

[p,gof] = fit(crosstimes(1,:)',crosstimes(2,:)','poly2');
% fitted = polyval(p,(0:1e-8:1e-5));

% figure();
% subplot(1,2,1);
% plot(tt_avg);
% subplot(1,2,2);
% hold on;
% scatter(crosstimes(1,:),crosstimes(2,:));
% plot(crosstimes(1,:),p(crosstimes(1,:)));

refgrid = (crosstimes(1,1):(tmeas(2)-tmeas(1)):crosstimes(1,end));
ygrid=(0:crosstimes(2,end)/(length(refgrid)-1):crosstimes(2,end));
xgrid = (-p.p2 + sqrt(p.p2^2-4*p.p1*(p.p3-ygrid)))/(2*p.p1);


res.xgrid = xgrid;
% save('resampling_points_6250pts_10us.mat','-struct','res');
% refsample = refgrid(2:end)-refgrid(1:end-1);
% xsample = xgrid(2:end)-xgrid(1:end-1);

% plot(refsample)
% hold on
% plot(xsample)

resampled = zeros(1,length(xgrid));
point = zeros(1,length(xgrid));
resampled_data = zeros(lenx,length(xgrid));
anchor = 210;

for i=1:length(xgrid)
    pointer=1;
    flag=0;
    while flag == 0
        if (refgrid(pointer) > xgrid(i) || pointer >= length(refgrid))
            flag = 1;
        else
            pointer=pointer+1;
        end        
    end
    
    resampled(i) = refgrid(pointer);
    point(i)=pointer;
    for k=1:lenx
        if pointer ==1
            resampled_data(k,i) = wmeas(k,1+anchor);
        elseif pointer == length(refgrid)
            resampled_data(k,i) = wmeas(j,anchor+length(refgrid));
        else
            resampled_data(k,i) = wmeas(k,pointer+anchor);
        end
    end
%     diff(i) = pointer
end

figure(1)
hold on;
plot(tmeas_resample,wmeas(1,anchor+1:anchor+1+length(xgrid)-1));
plot(tmeas_resample,resampled_data(1,:));
legend('Original','Resampled');
legend boxoff;
xlabel('Time (s)');
ylabel('Voltage (V)');
axis([0 0.9e-5 -2 2]);
% figure(2)

% realcrossing = zeros(1,cunt);


%%
[lenx,leny] = size(resampled_data);
tt_resample = zeros(lenx,57);
tmeas_resample = (tmeas(2)-tmeas(1))*(0:1:leny-1);
for j=1:lenx
    wmt_resample = resampled_data(j,:);
    crossingtime_resample = zeros(1,100);
    temp_resample = zeros(1,leny);
    count_resample = 1;

    for i=1:leny-1
        if wmt_resample(i) * wmt_resample(i+1) < 0
    %         if count ~= 0
    %         if tmeas(i) - crossingtime(count) > 0 
            count_resample = count_resample + 1;
            crossingtime_resample(count_resample) = tmeas_resample(i);
            temp_resample(i) = 1;
    %         end
    %         end
        end
    end

    crossingtime_resample = crossingtime_resample(2:count_resample);
    periodsraw_resample = crossingtime_resample(2:end)-crossingtime_resample(1:end-1);

    tempt_resample = [];
    for i=1:length(periodsraw_resample)
        if periodsraw_resample(i) > 0.5e-7
            tempt_resample = [tempt_resample periodsraw_resample(i)];
        end
    end

    tt_resample(j,:) = tempt_resample(2:58);
%     lent = length(tt);
end

tt_avg_resample = mean(tt_resample,1);

crosstimes_resample = zeros(2,length(tt_avg_resample));

for i=1:length(crosstimes_resample)
    crosstimes_resample(1,i) = sum(tt_avg_resample(1:i));
    crosstimes_resample(2,i) = i-1;
end

figure(1)
hold on
[per,F]=periodogram(wmeas(1,anchor+1:end),window(@rectwin,length(wmeas(1,anchor+1:end))),length(wmeas(1,anchor+1:end)),1/(tmeas(2)-tmeas(1)));
plot(F,10*log10(per))
[per_resample,F]=periodogram(resampled_data(1,:),window(@rectwin,length(resampled_data(1,:))),length(resampled_data(1,:)),1/(tmeas_resample(2)-tmeas_resample(1)));
plot(F,10*log10(per_resample))
axis([0 10e6 -120 -50])
legend('original','resample')
legend boxoff
% count = 1;

% for i=1:length(crossingtime)-1
%     if crossingtime(i+1)-crossingtime(i) > 100*(tmeas(2)-tmeas(1))
%         realcrossing(count) = crossingtime(i+1);
%         count = count + 1;
%     end
% end
% 
% realcrossing = realcrossing(1:count-1);