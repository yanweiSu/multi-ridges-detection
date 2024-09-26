clear;% close all;

load('./sampleTFR1.mat');
% load('./sampleTFR2.mat');

fs = basicTF.fs;
fr = basicTF.fr;
cALL = [];

%%
% Partition the TFR in time axis
tt = 0:fs:size(tfrsq,2);
% The first part's searching band (wider)
lowFreqINIT = 0.3; highFreqINIT = 4.0;

% Construct curves segment-by-segment
flag = 0;
cEND = [0; 0; 0];
for t = 1:length(tt)-1
    disp([num2str(t),'/', num2str(length(tt)),': ',num2str(lowFreqINIT), 'Hz - ', num2str(highFreqINIT),'Hz'])
    TFR = abs(tfrsq(:, tt(t)+1:tt(t+1)));
    tic

    [c0,c1,c2] = CurveMultiExt_init(TFR.', tfrsqtic, ...
        1.0, 0.8, 0.6, 800*fs, 600*fs, ...      % smooth(lambda) and similarity(mu) penalties
        2, 3, ...                               % multiples are 2 and 3
        lowFreqINIT/fs, highFreqINIT/fs, ...    % The fundamental's searching band
        round(0.2/fr), round(0.4/fr), ...       % The error(searching range) of the multiples
        flag, cEND);                            % The starting point 'cEND' of the extraction

    % The starting point for the next segment
    cEND = [c0(end); c1(end); c2(end)];
    % The searching band (finer)
    lowFreqINIT = tfrsqtic(c0(end))*fs - 0.2;
    highFreqINIT = tfrsqtic(c0(end))*fs + 1.0;

    cALL = [cALL; [c0 c1 c2]];  
    flag = 1;
    toc
end

%%
figure;
set(gcf,'Position',[100 50 1000 700]);
imageSQ((0:size(tfrsq,2)-1)./fs, tfrsqtic*fs, abs(tfrsq), 0.99);
axis xy; colormap(1-gray); colorbar
xlabel('time(sec)','FontSize',20);
ylabel('frequency(Hz)','FontSize',20);
ax = gca;
ax.FontSize = 20;
for k = 1:3
hold on;
plot((0:size(cALL,1)-1)./fs, tfrsqtic(cALL(:,k))*fs, 'r', 'LineWidth', 1.0)
end
ylim([0 20])
hold off;