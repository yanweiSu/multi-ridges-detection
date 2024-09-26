clear; close all;

% load('./sampleTFR1.mat');
load('./sampleTFR2.mat');

%% MultiCurve generic(2)
fs = basicTF.fs;
fr = basicTF.fr;
cALL = [];

% Partition the TFR in time axis
tt = 0:fs:size(tfrsq,2);
% The first part's searching band (wider)
lowFreqINIT = 0.3; highFreqINIT = 4.0;

% Construct curves segment-by-segment
flag = 0;
cEND = [0; 0];
for t = 1:length(tt)-1
    disp([num2str(t),'/', num2str(length(tt)),': ',num2str(lowFreqINIT), 'Hz - ', num2str(highFreqINIT),'Hz'])
    TFR = abs(tfrsq(:, tt(t)+1:tt(t+1)));
    
    tic
    [c0,c1] = CurveMultiExt_init_2curves(TFR.', tfrsqtic, ...
        1.0, 0.8, 800*fs, ...                   % smooth(lambda) and similarity(mu) penalties
        3, ...                                  % multiple is 3
        lowFreqINIT/fs, highFreqINIT/fs, ...    % The fundamental's searching band
        round(0.2/fr), ...                      % The error(searching range) of the multiples
        flag, cEND, 1.0);                       % The starting point 'cEND' of the extraction

    % The next starting point
    cEND = [c0(end); c1(end)];
    lowFreqINIT = tfrsqtic(c0(end))*fs - 1.0;
    highFreqINIT = tfrsqtic(c0(end))*fs + 1.0;
    cALL = [cALL; [c0 c1]];  
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
for k = 1:2
hold on;
plot((0:size(cALL,1)-1)./fs, tfrsqtic(cALL(:,k))*fs, 'r', 'LineWidth', 1.0)
end
ylim([0 20])
hold off;