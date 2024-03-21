function status = getStatus(targetDirectory)
    % Returns an structure with status of each worker.

    % Read all status files. Count occurances.
    statFiles = dir(strcat(targetDirectory, '/*.status'));
    % nstatus=cell(length(statFiles),1);
    status = struct;
    status.jobs = struct;
    status.nCompleted = 0;
    status.nFailed = 0;
    status.nCount = length(statFiles);

    for f = 1:length(statFiles)
        nstatus = fileread(fullfile(statFiles(f).folder, statFiles(f).name));

        if ~isfield(status.jobs, nstatus)
            status.jobs.(nstatus) = [f];
        else
            status.jobs.(nstatus) = [status.jobs.(nstatus), f];
        end

        if nstatus == "COMPLETED"; status.nCompleted = status.nCompleted +1; end
        if nstatus == "FAILED"; status.nFailed = status.nFailed + 1; end
    end

end
