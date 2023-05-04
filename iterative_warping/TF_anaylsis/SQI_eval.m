function [TakeorNot, entropy] = SQI_eval(x, fs ,len, pad)
    % cut the signal to segment of length 'len' with padding length 'pad'
    % (pad forward)
    
    threshold = 0.062;
    
    tmp  = buffer(x,len*fs,(len-pad)*fs)';
    tmp = tmp(len/pad:end,:);
    entropy_list = zeros([size(tmp,1),1]);
    TakeorNot = 1;
    for i = 1:size(tmp,1)
        [counts,edges] = histcounts(tmp(i,:), 'BinWidth', 1e3);
        counts = counts/size(tmp,2);
        entropy_list(i) = mean(-sum((counts+ 1e-10).*log(counts+ 1e-10),1)); % prevent Inf
    end
    if min(entropy_list) < threshold
        TakeorNot = 0;
    end
    entropy = min(entropy_list);
end