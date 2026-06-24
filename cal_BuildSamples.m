function smartSample = cal_BuildSamples(time)
% Build smart numerical integration samples
logOffsets = [0, 1e-7, 1e-3, 3e-3, 7e-3, ...
              1.5e-2, 3e-2, 5e-2, 1.5e-1, 0.45, 1.35];
smartSample = [];
for i = 1:length(time)-1
    dt = time(i+1) - time(i);
    idx = logOffsets < dt;
    smartSample = [smartSample; time(i) + logOffsets(idx)']; %#ok<AGROW>
end

smartSample = [smartSample; time(end)];
smartSample = unique(smartSample);
smartSample = sort(smartSample(:));
end
