%% written by Taehwan Kim (taehwan@berkeley.edu)
header;
c = 3e8;
% alpha = 10e6/(352/3e8);

averaged_per = zeros(5,3906);
averaged_solution = zeros(1,5);

distances=[10 11];

for j=1:2
    datapathname = strcat('~/isg/MOABB/lidar_experiments/experiments/data/laser/linewidth/cleo_dec9/result_',string(distances(j)*30),'m_100u_191mA_2.mat');
    sample = load(datapathname);

    wmeas = transpose(sample.wmeas(:,:));
    tmeas = transpose(sample.tmeas(2:end,:));

%     wstd = max(wmeas,[],2);
%     nor = (1/wstd(1))*eye(length(wstd));
%     wmeas = nor*wmeas;
    leng = length(tmeas)-1;
%     wmeas = wmeas(:,1:leng);
%     tmeas = tmeas(:,1:leng);
    wmeas = wmeas(:,end-leng+1:end);
    tmeas = tmeas(:,end-leng+1:end);
    
    [lenx, leny] = size(wmeas);

    [per,F]=periodogram(wmeas(1,:),window(@rectwin,length(wmeas(1,:))),length(wmeas(1,:)),1/(tmeas(2)-tmeas(1)));

    per_data = zeros(lenx,length(per));

    solutions_per = zeros(lenx,1);

    for i=1:lenx

        [per_temp,f_temp] = periodogram(wmeas(i,:),window(@rectwin,length(wmeas(i,:))),length(wmeas(i,:)),1/(tmeas(2)-tmeas(1)));
%         per_temp = per_temp(50:end-50);
%         f_temp = f_temp(50:end-50);
        per_data(i,:) = per_temp;

        [Y, I] = max(per_data(i,:));

        solutions_per(i) = f_temp(I);
        
    end

% alpha = mean(solutions_per)/((6.6/6.3)*352*10/8/3e8);

    averaged_per(j,:) = mean(per_data,1);
    averaged_solution(j) = mean(solutions_per);
%     pers = zeros(9,size(averaged_per,2));
%     pers(9,:) = averaged_per;

end
    
figure(1);
set(gcf,'DefaultAxesColorOrder',colororder);
hold on;
for i=1:2
    plot(f_temp,10*log10(averaged_per(i,:)));
end
% xlim([0 2e7]);

% p = polyfit(1.49*30*distances,averaged_solution,1);
% alpha=p(1)*3e8;
% 
% fbw=alpha*(tmeas(end)-tmeas(1));
% distance_firstpoint = 3e8*averaged_solution(1)/alpha;
% distances = 

% alpha*distance/3e8=frequency

% plot(f_temp,10*log10(pnoise_model_data(solution_lse(end,:),f_temp)));
% plot(f_temp,10*log10(averaged_per));
% plot(f_roi,(pnoise_model_data(solution_lse(1,:),f_roi)));
% plot(f_roi,(pnoise_model_data(solution_lse(10,:),f_roi)));
% plot(f_roi,(pnoise_model_data(solution_lse(20,:),f_roi)));
% std(solutions_per)/std(solutions_lse)
% xlim([0 2e7]);

% mean(solutions_per)*c/alpha/2
% mean(solutions_lse)*c/alpha/2
% 
% d_lse = solutions_lse*c/alpha/2;
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

