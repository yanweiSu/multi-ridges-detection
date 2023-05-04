function [sig, tfr, tfrtic, phi_value] = iterWarping(sig, basicTF, phi, I)
% This program iteratively warps on signal's fundamental, 1st harmonics, 2nd
% harmonics, ... and so on. The number of iterations is set by [I].
% [sig] is the signal
% [phi] is the pre-extracted fundamental phase
% [basicTF] is as following:


fs = basicTF.fs;
win = basicTF.win;
fr = basicTF.fr;
HighFreq = basicTF.HighFreq;
LowFreq = basicTF.LowFreq;
hop = 1;

% Amplitude (for reconstruction)
[h, ~, ~] = hermf(win,1,5);
h0 = h(floor(size(h,2)/2)+1);

phi_value = cell(I,1);
for i = 1:I
    disp(['warping ',num2str(i)]);
    %% construct "inverse of phi(psiInv)"
    M = ceil(phi(end)-phi(1))*fs;  % sampling on the fundamental phase
    % tau = 1/fs;
    val = linspace(phi(1), phi(end), M);
    psiInv = zeros(1,floor(M/i));
%     psiInv = zeros(1,floor(M/1));
    for k = 1:length(psiInv)
        [~, psiInv(k)] = min(abs(phi-val(i*k)));
%         [~, psiInv(k)] = min(abs(phi-val(k)));
    end

    %% Prepare for unwarping back (i.e. phi_hat in the note)
    phi_val = zeros(size(phi));
    val2 = val(i:i:i*floor(M/i));
%     val2 = val(1:1:floor(M));

    for k = 1:length(phi_val)
        [~, phi_val(k)] = min(abs(val2-phi(k))); % phi's value index
    end
    phi_value{i} = phi_val;

    %% Warping and SST on the warped signal
    sig = sig(psiInv);
    sig = sig - mean(sig);
    [~, ~, tfr, ~, tfrtic] = ConceFT_sqSTFT_C(sig, LowFreq, HighFreq, fr/fs, hop, win, 1, 5, 1, 1, 0);

    % TFR plot (For checking)
%     if i <= I
%         figure;
%         set(gcf,'Position',[100 50 1000 700]);
%         imageSQ((0:size(tfr,2)-1)./fs, tfrtic*fs, abs(tfr), 0.99);
%         axis xy; colormap(1-gray); colorbar
%         xlabel('time(sec)','FontSize',20);
%         ylabel('frequency(Hz)','FontSize',20);
%         ax = gca;
%         ax.FontSize = 20;
%     end

    %% extract the fundamental phase phi of the warped signal
    % The range of fundamental IF. Need to modify
    idx0 = find(tfrtic*fs>((i+1)-0.4) & tfrtic*fs<((i+1)+0.4));
%     idx0 = find(tfrtic*fs>((0+1)-0.3) & tfrtic*fs<((0+1)+0.3));
    [fund] = CurveExt(abs(tfr(idx0,:))', 0.6);
    fund = fund + idx0(1) - 1;
    tmp = Recon_sqSTFT_v2(tfr, tfrtic, fs/hop, fund, 0.5, h0);
    phi = unwrap(angle(tmp))/2/pi;  % This is the phase

end