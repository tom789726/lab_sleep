clc; clear; close all;

%% Loading data
% NPress data
filename = 'lab_sleep/{1AA7A696-78C6-45CD-B2FB-F9179E3CFA73}.SLP.csv'; 
nPress = dlmread(filename,','); 
pt = size(nPress,1);
fsamp = 25;

% Sleep event
filename = 'lab_sleep/{1AA7A696-78C6-45CD-B2FB-F9179E3CFA73}.SLP_OSA.csv'; 
E = dlmread(filename,','); 
Estart = E(:,1);
Edur = E(:,2);

%% Normalize to -1~1
nPress = nPress/max(abs(nPress));
time = 0:(1/fsamp):(1/fsamp)*(pt-1);

%% Crop data
bound = 0.0035;
nPress_pos = nPress.*(nPress>=bound);
nPress_neg = nPress.*(nPress<=(bound*-1));

%% Display

% figure
% plot(time,nPress); xlabel('Time(sec)'); grid on
% ax = gca; 
% ax.XAxis.Exponent = 0;
% % Mark apnea events
% hold on 
% y = ylim;
% plot([Estart Estart],[y(1) y(2)],'r-.');
% plot([Estart+Edur Estart+Edur],[y(1) y(2)],'b-.');
% 
% % obserce specific region
% xlim([12600 13600])


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
Vt = nPress_pos;
% Find zero crossing (count no. of cycles)

% % Testing
sig = nPress_pos(1:windowSz);
idx_zeros = find(sign(sig)==0);

sig_binary = sign(sig);
sig_pattern = sig_binary(1:end-1)*2+sig_binary(2:end);
idx_zCross = cat(1,[find(sig_pattern==1);find(sig_pattern==2)+1]); 
% See comm. sys.
% 00 = 0, 01 = 1
% 10 = 2, 11 = 3
idx_zCross = sort(idx_zCross);

figure
plot(time(1:windowSz),sig); xlabel('time(sec)');
hold on
plot(time(1+idx_zCross-1),sig(idx_zCross),'r*');

ylim([-0.5 3]);

% % END

%% Method 2: Envelope
% hold on
% nPress_env = zeros(size(nPress));
% [yUp,yLow] = envelope(nPress,fsamp,'peak');
% nPress_env = yUp;
% plot(time,yUp+0.1,'k-');

%% Display2

figure
plot(time,nPress_pos); xlabel('Time(sec)'); grid on
ax = gca; 
ax.XAxis.Exponent = 0;
% Mark apnea events
hold on 
y = ylim;
plot([Estart Estart],[y(1) y(2)],'r-.');
plot([Estart+Edur Estart+Edur],[y(1) y(2)],'b-.');

% obserce specific region
xlim([12600 13600])
% xlim([1 time(windowSz)]);
ylim([-0.5 1])
