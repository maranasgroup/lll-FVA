function rxnLink = getRxnLink(model, N, conComp, efmLoops)

if nargin < 3 || isempty(conComp)
    if size(N, 1) == sum(sum(model.S ~= 0, 1) > 1)
        conComp = -ones(size(model.S, 2), 1);
        conComp(sum(model.S ~= 0, 1) > 1) = connectVars(N);
    elseif size(N, 1) == size(model.S, 2)
        conComp = connectVars(N);
    else
        error('Nullspace matrix N has an incompatible size with model.S')
    end
end

rxnLink = sparse(size(model.S, 2), size(model.S, 2));
if nargin < 4 || isempty(efmLoops)
    p = pwd;
    efmToolpath = strsplit(which('CalculateFluxModes'), filesep);
    efmToolpath = strjoin(efmToolpath(1: end - 1), filesep);
    
    cd(efmToolpath)
    
    % [row, col, entry] = deal([]);
    % nEFM = 0;
    for jC = 1:max(conComp)
        try
            S = model.S(:, conComp == jC);
            S = S(any(S, 2), :);
            efms = CalculateFluxModes(full(S), double(model.lb(conComp == jC)<0), ...
                CreateFluxModeOpts('sign-only', true, 'level', 'WARNING'));
            pause(1e-4)
            efms = efms.efms;
            rxnJC = find(conComp == jC);
            for j = 1:numel(rxnJC)
                rxnLink(rxnJC(j), rxnJC) = any(efms(:, efms(j, :) ~= 0), 2)';
            end
        catch msg
            fprintf('Error encountered during calculation of EFMs:\n%s', getReport(msg))    
        end
    end
    cd(p)
else
    for j = 1:size(model.S, 2)
        rxnLink(j, :) = any(efmLoops(:, efmLoops(j, :) ~= 0), 2)';
    end
end


