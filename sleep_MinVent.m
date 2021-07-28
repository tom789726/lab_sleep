clc; clear; close all;

%% Loading data
filename = 'lab_sleep/1.csv'; % Change file name
NPress = dlmread(filename,',',[1 0 540000 0]); % Crop first row (label)
pt = size(NPress,1);

%% Normalize to -1~1
norm = NPress/max(abs(NPress));

%% Set time & Cut 5mins
time = 0:(1/25):(1/25)*(pt-1);
time = time';

startTime = 1500;
duration = 300;
% 5 mins
% startTime,duration: unit = sec
% pt,idx: unit = # (sample number)
idx = find((time-startTime)==0,1) : find((time - (startTime + duration) )==0,1); 
% find: return index of first non-zero element
% other feasible ways: sign(x): return 1 if x is positive

%% Version 1 (difference)
base = 0;
for i = 1:1500 % collect first 1 mins 
   base = base + abs(norm(i)); 
end

minVen_1 = zeros(pt,1);
minVen_1(1:1500) = base;
for i = 1500:pt
   minVen_1(i) = minVen_1(i-1) + abs(norm(i)-norm(i-1)) - abs(norm(i-1498)-norm(i-1499));
end

figure
subplot(2,1,1);
plot(time(idx),norm(idx)); xlabel('time(sec)');
title('NPress (Normalised)');
subplot(2,1,2);
plot(time(idx),minVen_1(idx)); xlabel('time(sec)');
legend('Ventilation (Sum of difference)');
title('Minute Ventilation (ver1)');

%% Version 2 (moving average)
window = 25*60; % 1 min
minVen_2 = movmean(norm,window);
minVen_2 = minVen_2/max(abs(minVen_2)); % normalise to -1~1


figure
subplot(3,1,1);
plot(time(idx),norm(idx)); xlabel('time(sec)');
title('NPress (Normalised)');
subplot(3,1,2);
plot(time(idx),minVen_2(idx)); xlabel('time(sec)');
title('Moving Average (60s,normalised)');
subplot(3,1,3);
title('Minute Ventilation (ver2)');
% Draw Histogram every second
height = 0;
for i = idx(1):25:idx(end-1)
    height = sum(minVen_2(i:i+25)) / 25; 
    % Note: absolute of sum =/= sum of absolutes
    height = abs(height);
    ax1 = area([time(i) time(i+25)],[height height],'FaceColor','r','EdgeColor','k'); hold on;
end
ax2 = plot(time(idx),norm(idx)+1); xlabel('time(sec)'); hold off
legend([ax1,ax2],{'Ventilation (Average in 1 sec)','NPress(norm)'});



%% Version 3 (volume & peak detection)
vol = cumtrapz(time,norm);
% figure
% findpeaks(vol(idx),'MinPeakProminence',0.015); title('Peak Detection(0.015)');
[pks,locs] = findpeaks(vol(idx),'MinPeakProminence',0.015);
% 0.5 -> 0.025 -> 0.015



figure
subplot(3,1,1);
plot(time(idx),norm(idx)); xlabel('time(sec)');
title('NPress (Normalised)');
subplot(3,1,2);
plot(time(idx),vol(idx)); xlabel('time(sec)');
title('Volume (integration)');
subplot(3,1,3);

% Draw Histogram by peak height & duration
height = 0;
locs_map = idx(locs); % Map from 0-7500 back to range of idx
idx_prev = idx(1);
for i = 1:size(locs_map,2)-1
    idx_next = locs_map(i);
    height = sum(vol(idx_prev:idx_next)) / (idx_next-idx_prev); % Divide by sample, not time
    height = abs(height);
    % Again, beware that absolute of sum =/= sum of absolutes
    ax1 = area([time(idx_prev) time(idx_next)],[height height],'FaceColor','r','EdgeColor','k');
    hold on;
    idx_prev = locs_map(i);
    idx_next = locs_map(i+1);
end
ax2 = plot(time(idx),norm(idx)+1); xlabel('time(sec)');
ax3 = plot(time(idx),vol(idx)); xlabel('time(sec)');
ax4 = plot(time(locs_map),vol(locs_map),'yo');
grid on; hold off;
legend([ax1,ax2,ax3,ax4],{'Minute Ventilation','NPress(Norm)','Volume','Peak'});
title('Minute Ventilation (ver3)');
