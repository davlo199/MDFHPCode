function status = waitFor(targetDir)
    % Wait for completion of jobs, return status object.

    % If waiting for more than this many seconds, die.
    timeoutTime = Inf;

    % How often to read status files.
    pollTime = 20;

    lastStatus = 0;
    startTime = tic;

    while 1
        status = getStatus(targetDir);
        % Make nice summary (might want to comment out for verboseness.
        if ~isequal(lastStatus, status)
            disp(status.jobs);
        end

        % Break loop if all done.
        if status.nCount == status.nCompleted + status.nFailed; break; end
        % Throw error if timeout.
        assert(toc(startTime) < timeoutTime, "Waiting for jobs timed out.");

        lastStatus = status;
        % Wait specified time.
        pause(pollTime);
    end

end
