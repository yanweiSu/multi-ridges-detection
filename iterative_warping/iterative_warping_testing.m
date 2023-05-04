clear;
close all;
addpath('.\TF_anaylsis');
% load('allsubLabel.mat');

%% Here is the input signal, assign it to be 'sig'
% load('allsubPPG.mat');
% load('allsubLabel.mat');
% sub = 8;
% PPG_sub8 = allsubPPG{sub};
load("normalPPG.mat");
% example:
ep = 648;   % N2     %43, 387, 513, 284
% ep = 549;   % REM
% ep = 290;   % N3
% ep = 63;    % AW
fs = 100;
len_epoch = 30; % sec
len_sig = 90; %sec

% Buffer
signal = buffer(PPG_sub8, len_sig*fs, (len_sig-len_epoch)*fs); % sampling rate is 100 Hz
signal = signal(:, len_sig/len_epoch:end); % The first one's labeling=2. The final one is num_epochs-1.
slabel = ceil(len_sig/len_epoch/2); % The first one's label
% disp(allsubLabel{sub}(ep));

% IFcells = cell(length(allsubLabel{sub}),1); % store the IF curves
% features_PHI = cell(length(allsubLabel{sub}),1); % store PHIcurves

% parfor ep = slabel : length(allsubLabel{sub})-(slabel-1)
dws = 1;
sig = signal(:,ep-slabel+1);
sig = sig - mean(sig);
% sig = resample(sig, fs/dws, fs);
fs = fs/dws;


% Quality
% qual = ppgSQI(sig(len_epoch*fs+1:2*len_epoch*fs), fs);
% disp(qual);

%% First get the fundamental component: TF pararmeters setting and prepare the TFR
fr = 0.01;
win = fs*8+1;
hop = 1;
HighFreq = 10/fs;
LowFreq = 0.1/fs;
[~, ~, tfrsq, ~, tfrsqtic] = ConceFT_sqSTFT_C(sig, LowFreq, HighFreq, fr/fs, hop, win, 1, 5, 1, 1, 0);

% Amplitude (for reconstruction)
[h, ~, ~] = hermf(win,1,5);
h0 = h(floor(size(h,2)/2)+1);

figure; plot((0:length(sig)-1)./fs, sig, 'LineWidth', 1.2); title('signal');
ax = gca; ax.FontSize = 16; set(gca, 'XTickLabel', []); set(gca, 'YTickLabel', []);
hold on;
% [b, a] = butter(5, [0.5, 12] / (fs / 2));
% sig = filtfilt(b, a, sig);
% plot((0:length(sig)-1)./fs, sig, 'LineWidth', 1.2);
ylim([-3 3]*1e4)
hold off;

% TFR plot
figure; set(gcf,'Position',[100 50 1000 700]);
subplot(1,2,1);
imageSQ((0:size(tfrsq,2)-1)./fs, tfrsqtic*fs, abs(tfrsq), 0.99);
axis xy; colormap(1-gray);% colorbar
xlabel('time(sec)','FontSize',20); ylabel('frequency(Hz)','FontSize',20);
ax = gca; ax.FontSize = 20;
title(['2ndSST: window = ',num2str((win-1)/fs),' sec']);

% Get fundamental's phase [phi_fund]
idx0 = find(tfrsqtic*fs>0.5 & tfrsqtic*fs<1.9);
[fund] = CurveExt(abs(tfrsq(idx0,:))', 0.8);
fund = fund + idx0(1) - 1;
hold on;
plot((0:length(fund)-1)/fs, fs*tfrsqtic(fund), 'r');
hold off;
tmp = Recon_sqSTFT_v2(tfrsq, tfrsqtic, fs/hop, fund, 0.5, h0);
phi_fund = unwrap(angle(tmp))/2/pi;  % This is the phase of the fundamental

% subplot(2,1,2);
% plot((1:length(phi_fund)-1)./fs, diff(phi_fund)*fs);
% % phi_fund = medfilt1(medfilt1(phi_fund, 100), 100);
% % plot(diff(phi_fund)*fs); hold off;
% ylim([0.5 1.5])
% xlabel('time(sec)','FontSize',20); ylabel('frequency(Hz)','FontSize',20); title('\phi_1^\prime');
% ax = gca; ax.FontSize = 20;

%% Iterative warping
I = 3; % the number of iterations
basicTF.fs = fs;
basicTF.fr = fr;
basicTF.win = fs*7+1;
basicTF.HighFreq = 10.0/fs;
basicTF.LowFreq = 0.1/fs;

[sig_warped, tfr_warped, tfrtic_warped, phi_value] = iterWarping(sig, basicTF, phi_fund, I);

% Amplitude (for reconstruction)
[h, ~, ~] = hermf(basicTF.win,1,5);
h0_warp = h(floor(size(h,2)/2)+1);

%% CurveExt on the warped TFR (Use fundamental-informed version)
numH = 6;
% % 1. Get the fundamental
% idx0 = find(tfrtic_warped*fs>0.5 & tfrtic_warped*fs<1.9);
% [c] = CurveExt(abs(tfr_warped(idx0,:))', 3.0);
% c = c + idx0(1) - 1;
% 
% cALL = zeros(length(c), 7);
% cALL(:,1) = c;
% 
% % 2. Multiple curve extraction
% tic;
% [cALL(:,2), cALL(:,3), cALL(:,4)] = CurveExt_multi2_ver0(abs(tfr_warped).', fs*tfrtic_warped, c, ...
%     3, 3, 3, ... %lambda_k: smooth penalty
%     10, 8, 6, ... %mu_k: similarity penalty
%     2, 3, 4, ... %2nd, 3rd, 4th harmonics
%     5, 5, 5);   %bandwidth = \pm5*fr = \pm0.1Hz
% toc;
% 
% tic;
% [cALL(:,5), cALL(:,6), cALL(:,7)] = CurveExt_multi2_ver0(abs(tfr_warped).', fs*tfrtic_warped, c, ...
%     3, 3, 3, ...
%     5, 4, 3, ...
%     5, 6, 7, ...
%     6, 6, 6);
% toc;
% 
% % Plot for checking
% figure;
% set(gcf,'Position',[100 50 1000 700]);
% imageSQ((0:size(tfr_warped,2)-1)./fs, tfrtic_warped*fs, abs(tfr_warped), 0.99);
% axis xy; colormap(1-gray); colorbar
% xlabel('time(sec)','FontSize',20);
% ylabel('frequency(Hz)','FontSize',20);
% ax = gca;
% ax.FontSize = 20;
% for k = 1:7
%     hold on;
%     plot((0:length(c)-1)./fs, fs*tfrtic_warped(cALL(:,k)), 'r');
%     hold off;
% end

% %% CurveExt on the warped TFR: MultiCurve_generic(3)
% fr = basicTF.fr;
% tic
% cALL = []; cEND = [1; 1; 1];
% tt = [0:0.2*fs:size(tfr_warped,2) size(tfr_warped,2)];
% lowFreqINIT = 0.7; highFreqINIT = 1.3; flag = 0;
% for t = 1:length(tt)-1
%     TFR = abs(tfr_warped(:, tt(t)+1:tt(t+1)));
%     [c0,c1,c2] = CurveMultiExt_init(TFR.', tfrtic_warped, ...
%         1.0, 0.8, 0.6, 1500*fs, 1400*fs, ...
%         2, 3, ...
%         lowFreqINIT/fs, highFreqINIT/fs, round(0.2/fr), round(0.3/fr), ...
%         flag, cEND, 1.0);
%     disp([num2str(t),'/', num2str(length(tt)),': ',num2str(lowFreqINIT), 'Hz - ', num2str(highFreqINIT),'Hz'])
%     cEND = [c0(end); c1(end); c2(end)];
%     lowFreqINIT = tfrtic_warped(c0(end))*fs - 0.2;
%     highFreqINIT = tfrtic_warped(c0(end))*fs + 0.2;
%     cALL = [cALL; [c0 c1 c2]];  
%     flag = 1;
% end
% toc;
% lowFreqINIT = tfrtic_warped(cALL(1,1))*fs - 0.3;
% highFreqINIT = tfrtic_warped(cALL(1,1))*fs + 0.3;
% 
% for h = 4:6
% [cALL(:,h)] = CurveMultiExt_withFund(abs(tfr_warped).', tfrtic_warped, ...
%     cALL(:,1), 0.6, 1000*fs, h, round(0.4/basicTF.fr));
% end
% % delta = 11+(1:8).';
% % harmonics = IFext_multi(tfr_warped, tfrtic_warped, fs, 1, 6, delta); % This is the extract one-by-one version
% toc
% harmonics = cALL;

% CurveExt
cALL = CurveExt(abs(tfr_warped).', 1.0);
for h = 2:numH
    fprintf("%d\n", h);
    cALL = [cALL CurveMultiExt_withFund(abs(tfr_warped).', tfrtic_warped, ...
        cALL(:,1), 0.6, 0*fs, h, round(0.4/basicTF.fr))];
    % Usage of CurveMultiExt_withFund: See the C file
end

% Plot for checking
subplot(1,2,2); set(gcf,'Position',[100 50 1000 700]);
imageSQ((0:size(tfr_warped,2)-1)./fs, tfrtic_warped*fs, abs(tfr_warped), 0.99);
axis xy; colormap(1-gray); %colorbar
xlabel('time(sec)','FontSize',20); ylabel('frequency(Hz)','FontSize',20);
ax = gca; ax.FontSize = 20;
for k = 1:numH
    hold on;
    % plot((tt(1)+1:tt(end))./fs, fs*tfrtic_warped(cALL(:,k)), 'r');
    plot((0:length(sig_warped)-1)./fs, fs*tfrtic_warped(cALL(:,k)), 'r');
    hold off;
end
title('Warped TFR');

%% Recon by SAMD and superposition
fundAM = []; fundFM = []; recon = []; super = [];
for h = 1:numH
    tmp = Recon_sqSTFT_v2(tfr_warped, tfrtic_warped, fs, cALL(:,h), 0.5, h0_warp);
    fundFM = [fundFM; unwrap(angle(tmp))];
    fundAM = [fundAM; abs(tmp)];
    recon = [recon; tmp];
    if h==1; super = tmp; else; super = super + tmp; end
end
super = real(super);

[fest, f1est, coeff] = SAMD_nophase(sig_warped.', fundAM(1,:), fundFM(1,:), numH, 2, fundFM);
% [fest, f1est, coeff] = SAMD(sig_warped.', fundAM, fundFM(1,:), 6, 2);
% figure; plot(sig_warped); hold on; plot(f1est); hold off;
% figure; plot(sig); hold on; plot(fest); hold off;

for i = I:-1:1 % Unwarp back
    fest = fest(phi_value{i});
    super = super(phi_value{i});
    fundFM = fundFM(:, phi_value{i});
    recon = recon(:, phi_value{i});
end

tt = 80*fs-5*fs+1:80*fs+5*fs; % tt = (fs*1+1):(length(sig)-fs*1);
disp(['SAMD_Error = ', num2str(norm(fest(tt).'-sig(tt))/norm(sig(tt)))]);
disp(['superpose_Error = ', num2str(norm(super(tt).'-sig(tt))/norm(sig(tt)))]);

figure; plot((0:length(sig)-1)./fs, sig); hold on; plot((0:length(sig)-1)./fs, super); hold off;
ax = gca; ax.FontSize = 16;
% set(gca, 'XTickLabel', []); set(gca, 'YTickLabel', []);
ylim([-3 3]*1e4)
title('Recon')

% figure;
% subplot(4,1,1); plot(tt./fs, sig(tt));
% subplot(4,1,3); plot(tt./fs, mod(fundFM(2,tt),2*pi));
% subplot(4,1,4); plot(tt./fs, (2*fundFM(1,tt)-fundFM(2,tt))./2/pi);% ylim([-0.5 0])
% subplot(4,1,2); plot(tt./fs, mod(fundFM(1,tt),2*pi));
% % hold on; plot(tt(1:end-1)./fs, (fs.*diff(fundFM(1,tt)))/max(fs.*diff(fundFM(1,tt)))*max(2*pi)); hold off;

% FM = mod(fundFM(1:6,:), 2*pi);
% rng = range(FM(1:6,:));
% 
% figure;
% plot(tt./fs, (sig(tt)-min(sig(tt)))./max(sig(tt)-min(sig(tt)))*2*pi, 'k'); hold on;
% ind = [tt(1)];
% peaks = [];
% for t = tt(1:end-1)
%     if FM(1,t)>4 && FM(1,t+1)<4
%         ind = [ind t];
%         [~,pks] = min(rng(ind(end-1):t));
%         peaks = [peaks pks+ind(end-1)-1];
%     end
% end
% plot(peaks./fs, (sig(peaks)-min(sig(tt)))./max(sig(tt)-min(sig(tt)))*2*pi, 'o');
% % plot(tt./fs, rng(tt));
% 
% for hh = 1:numH
% plot(tt./fs, mod(fundFM(hh,tt),2*pi),'.','MarkerSize',20); hold on;
% end
% hold off;
% % plot(tt./fs, mod(fundFM(2,tt),2*pi)); hold on;
% % plot(tt./fs, mod(fundFM(3,tt),2*pi)); hold off;
% 
% figure;
% for hh = 1:3
%     plot(tt(1:end-1)./fs, fs.*diff(fundFM(hh,tt))/2/pi);
%     hold on;
% end
% hold off;

%% PHI
PHIcurves = zeros(6, size(recon,2)); % minus 1 stands for the fundamental components
figure;
for j = 2:6
    if j >= 2
        PHI = ((recon(1,:)./abs(recon(1,:))).^j).*conj(recon(j,:)./abs(recon(j,:)));
%         PHI = mod((j*fundFM(1,:)-fundFM(j,:))./2/pi, 1);
    else
        PHI = recon(1,:);
    end
    PHI = unwrap(angle(PHI))/2/pi;
    if j == 1
        FM1 = PHI;
        X = [ones(length(FM1),1) (0:length(FM1)-1).'./fs];
        b = X\(FM1.');
        fitt = b(1)+b(2)*(0:length(FM1)-1)./fs;
        PHI = FM1-fitt;
    end

    % Phase fixing
%     PHI = fixing(PHI, fs);
%     if j==3
%         for tt = 1:length(PHI)
%             if ~isnan(PHI(tt)) && PHI(tt)<0.0
%                 PHI(tt) = PHI(tt) + 1.0;
%             end
%         end
%     end
    
    plot((0:length(PHI)-1)./fs, PHI, 'LineWidth', 1.2, 'LineStyle', '-'); hold on;
    PHIcurves(j,:) = PHI;
end
% plot((0:length(real(super))-1)./fs, real(super)/max(abs(super)), 'k', 'LineStyle', '-');
hold off; ylim([-1 1])
title('PHI curves');

%% AHI
figure;
AHIcurves = zeros(6,size(recon,2));
for j = 1:6
    AHIcurves(j,:) = abs(recon(j,:));
    if j==1
%     plot(tt./fs, ...
%         AHIcurves(j,tt), 'LineWidth', 1.2, 'LineStyle', '-'); hold on;
%         abs(recon(j,:)), 'LineWidth', 1.2, 'LineStyle', '-'); hold on;
    else
%     plot((0:size(AHIcurves,2)-1)./fs, ...
%         (AHIcurves(j,:)./AHIcurves(1,:)), ㄏ'LineWidth', 1.2, 'LineStyle', '-'); hold on;
    plot((0:length(sig)-1)./fs, ...
        (AHIcurves(j,:)./AHIcurves(1,:)), 'LineWidth', 1.2, 'LineStyle', '-'); hold on;
%         abs(recon(j,:)), 'LineWidth', 1.2, 'LineStyle', '-'); hold on;
    end
end
% plot((0:length(PHI)-1)./fs, PHIcurves(2,:), 'LineWidth', 1.2, 'LineStyle', '-');
% hold on;
% plot((0:length(sig)-1)./fs, real(super)/1.2/max(abs(super)), 'k', 'LineStyle', '-');
% hold off;
% % xline(peaks./fs);
% xlim([tt(1) tt(end)]./fs);