function results = collectResults(workDir)
    % Returns outputs re-glued together.
    numResults = length(dir(strcat(workDir, '/*output.mat')));
    results = [];

    % Load outputs variables.
    for i = 1:numResults
        results(:, i) = load(fullfile(workDir, strcat("worker", num2str(i), "output.mat")), 'returnValues').returnValues;
    end

end
