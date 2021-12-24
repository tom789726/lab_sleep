clc; clear; close all;

%% Parameters
% time_op = 10400;
% time_ed = 11400;
time_op = 10750;
time_ed = 11000;
Edur_min = 10; % sec

%% Loading data
% NPress data
filename = 'lab_sleep/{3A1A3189-1099-442A-80B7-F16AC41E69CC}.SLP.csv'; 
nPress = dlmread(filename,','); 
pt = size(nPress,1);
fsamp = 25;

% Sleep event
filename = 'lab_sleep/{3A1A3189-1099-442A-80B7-F16AC41E69CC}.SLP_OSA.csv'; 
E = dlmread(filename,','); 
Estart = E(:,1);
Edur = E(:,2);


%% Normalize to -1~1
nPress = nPress/max(abs(nPress));
time = 0:(1/fsamp):(1/fsamp)*(pt-1);

%% Transform pressure to flow
flow = sqrt(nPress);
bound_f = 0.04;
flow_pos = flow.*(flow>=bound_f);

%% Crop data (nPress)
bound = 0.0025;
nPress_pos = nPress.*(nPress>=bound);
nPress_neg = nPress.*(nPress<=(bound*-1));

%% Display (raw signal)

% figure
% plot(time,nPress); 
% % plot(time,flow);
% xlabel('Time(sec)'); ylabel('NPress'); grid on
% ax = gca; 
% ax.XAxis.Exponent = 0;
% % Mark apnea events
% hold on 
% y = ylim;
% plot([Estart Estart],[y(1) y(2)],'r-.');
% plot([Estart+Edur Estart+Edur],[y(1) y(2)],'b-.');
% 
% 
% % obserce specific region
% xlim([time_op time_ed])


%% Method 1: sliding window (60s)
windowSz = 60*fsamp;
stepSz = fsamp;

% nPress_slide = zeros(size(nPress));
% nPress_slide(1:windowSz) = sum(abs(nPress(1:windowSz)))/(2*windowSz);
% for i = 0:stepSz:pt
%     j = i+1;
%     if (j+windowSz <= pt)
%         nPress_slide(j:j+windowSz) = sum(abs(nPress(j:j+windowSz)))/(2*windowSz);
%     end
% end
% 
% hold on
% plot(time,nPress_slide+0.5,'g-');

%% 1.1: Tidal Volume

% Find zero crossing pts (+ count no. of cycles)
sig = nPress_pos;
% sig = flow_pos;
idx_zeros = find(sign(sig)==0);

sig_binary = sign(sig);
sig_pattern = sig_binary(1:end-1)*2+sig_binary(2:end);
idx_zCross = cat(1,[find(sig_pattern==1);find(sig_pattern==2)+1]); 
idx_zCross = sort(idx_zCross);
% See comm. sys.
% 00 = 0, 01 = 1
% 10 = 2, 11 = 3
disp("Total no. of cycles = " + (size(idx_zCross,1)-1));

% Find Tidal Volume
Vt = zeros(size(nPress_pos));
idx_vt = [idx_zCross(1:end-1),idx_zCross(2:end)];
% Vt_temp = double.empty();
for i = 1:size(idx_vt,1)
   Vt_dur = time(idx_vt(i,2)) - time(idx_vt(i,1));
   Vt(idx_vt(i,1):idx_vt(i,2)) = sum(sig(idx_vt(i,1):idx_vt(i,2))) / (Vt_dur/fsamp);
end

figure
ax1 = plot(time,sig); xlabel('time(sec)');
hold on
ax2 = plot(time(1+idx_zCross-1),sig(idx_zCross),'mx');
hold on
ax3 = plot(time,Vt/max(Vt(:)),'k-');

ax = gca; 
ax.XAxis.Exponent = 0;
title('Tidal Volume (Amount/Cycle)');


xlim([time_op time_ed])
% xlim([9160 9460])
% xlim([1 120])
ylim([-0.5 1])

hold on 
y = ylim;
plot([Estart Estart],[y(1) y(2)],'r-.');
plot([Estart+Edur Estart+Edur],[y(1) y(2)],'b-.');
legend([ax1,ax2,ax3],{'NPress(cropped)','zero crossings','Tidal Volume'});


%% 1.1.1 Tidal Volume with interpolation
% Edur_min: 判斷零點之間是否event

%% 1.1.2 Tidal Volume (Add negative signals)
% Find zero crossing pts (+ count no. of cycles)
sig_n = nPress_neg.*(-1);
% sig = flow_pos;
idx_zeros = find(sign(sig_n)==0);

sig_binary = sign(sig_n);
sig_pattern = sig_binary(1:end-1)*2+sig_binary(2:end);
idx_zCross = cat(1,[find(sig_pattern==1);find(sig_pattern==2)+1]); 
idx_zCross = sort(idx_zCross);
% See comm. sys.
% 00 = 0, 01 = 1
% 10 = 2, 11 = 3
disp("Total no. of cycles = " + (size(idx_zCross,1)-1));

% Find Tidal Volume

% **Do not clear memory of Vt to success the data from postive signals**
% Vt = zeros(size(nPress_neg));

idx_vt_n = [idx_zCross(1:end-1),idx_zCross(2:end)];
% Vt_temp = double.empty();
for i = 1:size(idx_vt_n,1)
   Vt_dur = time(idx_vt_n(i,2)) - time(idx_vt_n(i,1));
   Vt(idx_vt_n(i,1):idx_vt_n(i,2)) = Vt(idx_vt_n(i,1):idx_vt_n(i,2)) + sum(sig_n(idx_vt_n(i,1):idx_vt_n(i,2))) / (Vt_dur/fsamp);
end

figure
ax1 = plot(time,nPress); xlabel('time(sec)');
% hold on
% ax2 = plot(time(1+idx_zCross-1),nPress(idx_zCross),'mx');
hold on
ax3 = plot(time,Vt/max(Vt(:))+0.2,'k-');
% ax3 = plot(time,Vt/max(Vt(:)),'k-');

ax = gca; 
ax.XAxis.Exponent = 0;
title('Tidal Volume (Amount/Cycle)');


xlim([time_op time_ed])
% xlim([9160 9460])
% xlim([1 120])
ylim([-0.5 1])

hold on 
y = ylim;
plot([Estart Estart],[y(1) y(2)],'r-.');
plot([Estart+Edur Estart+Edur],[y(1) y(2)],'b-.');
% legend([ax1,ax2,ax3],{'NPress(cropped)','zero crossings','Tidal Volume'});
legend([ax1,ax3],{'NPress(raw)','Tidal Volume'});

%% 1.2 Tidal Volume * Respiratory Rate , simple sliding window (60sec)
sig = nPress_pos;
% sig = flow_pos;
Vt_window = zeros(size(sig));
idx_window = (1:windowSz);

RR_start = find(idx_zCross>=idx_window(1),1);
RR_end = find(idx_zCross>idx_window(end),1)-1;
RR = RR_end-RR_start;
disp("No. of respiratory cycles in this window = " + RR);

Vt_window(idx_window) = sum(Vt(idx_window))/60 * RR; % 前60秒(第一分鐘)用同一個數值

count = 0;
for i = windowSz:stepSz:pt 
    % Sliding window:
    % window size = 60sec (25Hz) , step size = 1sec (25Hz)
    
    idx_window = (1+stepSz*count:windowSz+stepSz*count);
    
    RR_start = find(idx_zCross>=idx_window(1),1);
    RR_end = find(idx_zCross>idx_window(end),1)-1;
    if (isempty(RR_end)) % Index Touches End of samples
       RR = idx_zCross(end) - idx_zCross(RR_start);
    else
        RR = RR_end-RR_start;
    end
    
    % disp("No. of respiratory cycles in this window = " + RR);
    
    % 1格 = 1秒 (一次更新25Hz資料點) , 由第60-61秒開始更新
    Vt_window(windowSz+stepSz*(count-1):i) = sum(Vt(idx_window))/60 * RR;
    
    count = count+1;
end



figure
ax1 = plot(time,sig); xlabel('time(sec)');
hold on
ax2 = plot(time,Vt_window/max(Vt_window),'k-');

ax = gca; 
ax.XAxis.Exponent = 0;
xlim([time_op time_ed])
% xlim([1 120])
ylim([-0.1 1])

hold on 
y = ylim;
plot([Estart Estart],[y(1) y(2)],'r-.');
plot([Estart+Edur Estart+Edur],[y(1) y(2)],'b-.');

title('Minute Ventilation');
legend([ax1,ax2],{'NPress(cropped)','Minute Ventilation'});
% legend([ax1,ax2],{'Flow(cropped)','Minute Ventilation'});

%% Display2 (cropped signal)

% figure
% plot(time,nPress_pos); xlabel('Time(sec)'); grid on
% ax = gca; 
% ax.XAxis.Exponent = 0;
% % Mark apnea events
% hold on 
% y = ylim;
% plot([Estart Estart],[y(1) y(2)],'r-.');
% plot([Estart+Edur Estart+Edur],[y(1) y(2)],'b-.');
% 
% % obserce specific region
% xlim([time_op time_ed])
% ylim([-0.2 0.5])
