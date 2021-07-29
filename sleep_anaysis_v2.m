clc; clear; close all;
%% Loading data
filename = 'data/DATA_sleep_analysis/ªL¯Î¹| 409_428 18736180 20190214.txt'; % Change file name
var = dlmread(filename,'\n',1,0);

%% Initialization
samplerate = 1/25;
var = var';
len = size(var,1);
time = 0:samplerate:(len-1)*samplerate;
time = time';

startTime = 0;
stopTime = 300; % adjusting interesting area
idx = find(sign(time-startTime)+1, 1) : find(sign(time-stopTime)+1, 1);
idx = idx';

%% Volume
int_var = cumtrapz(time,var);

figure
subplot(2,1,1)
plot(time(idx),var(idx),'k');
grid on;
xlabel('Time (s)');
ylabel('Flow');
xlim([min(time(idx)),max(time(idx))]);
subplot(2,1,2)
plot(time(idx),int_var(idx),'b');
grid on;
xlabel('Time (s)');
ylabel('Volume');
xlim([min(time(idx)),max(time(idx))]);

%% Minute ventilation
neg_var = -var;
[pks,locs] = findpeaks(neg_var(idx),'MinPeakProminence',0.0001);
% figure('Position',[1200 500 700 500],'Name','Peaks','NumberTitle','off')
figure
findpeaks(var(idx),time(idx),'MinPeakProminence',0.0001);hold on;
plot(time(locs),-pks,'ro');hold off

% section integral
mv_sec_var = zeros(size(locs));
figure('Position',[700 50 600 500])
for i = 1:size(locs,1)-1
    sec_var = var(locs(i):locs(i+1));
    int_sec_var = cumtrapz(time(1:size(sec_var),1),sec_var);
    mv_sec_var(i) = (max(sec_var)-sec_var(1))*60/(time(locs(i+1))-time(locs(i)));
    ax1 = area([time(locs(i)) time(locs(i)) time(locs(i+1)) time(locs(i+1))],[0 mv_sec_var(i) mv_sec_var(i) 0],'FaceColor',[0.8500 0.3250 0.0980]);hold on
%     plot(time(locs(i):locs(i+1)),sec_var,'b',time(locs(i):locs(i+1)),int_sec_var,'k');hold on
end
ax2 = plot(time(idx),var(idx),'b','LineWidth',2);
ax3 = plot(time(idx),int_var(idx),'r','LineWidth',2);
grid on;hold off
xlabel('Time(s)');
legend([ax1,ax2,ax3],{'min ventilation','flow','volume'});

%% Save file
[filename,filepath] = uiputfile({'*.csv'},'Àx¦sÀÉ®×',[filename(1:end-4),'.csv']);
idx = find(sign(time-startTime)+1, 1) : find(sign(time-stopTime)+1, 1);
str = strcat(filepath,filename);

data = [var(idx),int_var(idx)];
data2 = [time(locs),mv_sec_var];
header = array2table([data(1:size(data2,1),:),data2],'VariableNames',{'Flow','Volumn','Inital_location','Min_ventilation'});
writetable(header,str);
dlmwrite(str, data(size(data2,1):end,:),'-append');