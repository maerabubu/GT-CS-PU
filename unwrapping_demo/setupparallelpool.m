function setupparallelpool
%SETUPPARALLELPOOL Summary of this function goes here
%   Detailed explanation goes here
if isempty(gcp('nocreate'))
    nCpus = feature('numCores');
    parpool('local', nCpus);
end
end

