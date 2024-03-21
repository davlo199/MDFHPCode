% This script is designed to be called from slurmWrapper.sl
% 
% Loads broadcast and independant variables, then calls function and saves
% output.

% Get env variables.
workerID = getenv('SP_WORKER_ID');
cd(getenv('SP_WORKDIR'));

% Statusfile name.
statusFile = strcat('worker', workerID, '.status');

try
    path(pathdef);
    load('functionHandle.mat', 'functionHandle');
    load(strcat('worker', workerID, 'index.mat'), 'indexedVariable');
    fprintf(fopen(statusFile, 'w'), 'RUNNING');

    % Run
    returnValues = functionHandle(indexedVariable);
    save(strcat('worker', workerID, 'output.mat'), 'returnValues');
catch e
    fprintf(fopen(statusFile, 'w'), 'FAILED');
    error(e.message);
end

fprintf(fopen(statusFile, 'w'), 'COMPLETED');
