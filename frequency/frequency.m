%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Single-Tone Frequency Measurement          %
%              with MATLAB Implementation              %
%                                                      %
% Author: Ph.D. M.Sc. Eng. Hristo Zhivomirov  04/07/16 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f0 = frequency(x, fs)

% remove the DC component
x = x - mean(x);

% signal length
xlen = length(x);

% window generation
w = blackmanharris(xlen);

% find signal periodogram
Xm = periodogram(x, w, xlen)';

% coarse search for the signal frequency index via MLE
[~, indprim] = max(Xm);

% fine search for the signal frequency index via WAE
step = 4;
inds = (indprim-step:1:indprim+step);
ind = sum(Xm(inds).*inds)/sum(Xm(inds));

% signal frequency estimation
f0 = (ind-1)*fs/xlen;

end