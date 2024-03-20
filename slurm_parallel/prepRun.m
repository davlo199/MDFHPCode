function prepRun(workDir, functionHandle, indexedVariables)
    % prepRun

    % Creates .mat files required for matlab wrapper.

    % Save broadcastVariables and path
    save(fullfile(workDir, "functionHandle.mat"), 'functionHandle', '-mat');
    savepath(fullfile(workDir, 'pathdef.m'));

    % Save index/slice variables.
    for i = indexedVariables
        % Variables unique to this job
        indexedVariable = i;

        % Save for index;
        save(fullfile(workDir, strcat("worker", num2str(i), "index.mat")), 'indexedVariable', '-mat');
        fprintf(fopen(fullfile(workDir, strcat('worker', num2str(i), ".status")), 'w'), "PENDING");
    end
end
