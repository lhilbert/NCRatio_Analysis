function [threshold,metric] = otsuLimit(intValues)
% [threshold,metric] = otsuLimit(intValues)
%
% Finds Otsu threshold based on intensity values (intValues)
%
% Returns optimal threshold and Otsu metric for optimal threshold

intValues = intValues(:);

maxInt = max(intValues);
binEdges = 0:1:maxInt;
    
counts = histc(intValues,binEdges);
    
counts = counts(1:end-1);
binEdges = binEdges(1:end-1);
    
binEdges = (double(binEdges(counts>0))).';
counts = double(counts(counts>0));
    
counts = counts./max(counts);

numSupport = numel(binEdges);

pp = counts / sum(counts);
sigma_b = zeros(1,numSupport);

for tt = 1:numSupport
   q_L = sum(pp(1:tt));
   q_H = sum(pp(tt+1:end));
   mu_L = sum(pp(1:tt) .* binEdges(1:tt)) ./ q_L;
   mu_H = sum(pp(tt+1:end) .* binEdges(tt+1:numSupport)) ./ q_H;
   sigma_b(tt) = q_L.*q_H.*(mu_L-mu_H)^2;
end

[metric,maxInd] = max(sigma_b);

threshold = binEdges(maxInd);

end