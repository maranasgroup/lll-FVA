fprintf('Load the model for testing...\n')
p = which('iIT341.xml');
model = readCbModel(p);
% set bounds
ex = sum(model.S ~= 0, 1) <= 1;
model.lb(model.lb < 0 & ~ex') = -1000;
model.ub(model.ub > 0 & ~ex') = 1000;
model.lb(model.lb < -1000) = -1000;
model.ub(model.ub > 1000) = 1000;

methods = {'Original ll-FVA', 'll-FVA with Fast-SNP preprocessing', ...
    'Localized ll-FVA (nullspace-based)', 'Localized ll-FVA (EFM-based)'};
meth2test = 1:4;

fprintf('Start testing...\n')
rxnBiomass = 'BIOMASS_HP_published';
model.c(:) = 0;
model = changeObjective(model, rxnBiomass, 1);
for j = 1:numel(meth2test)
    t0 = tic;
    % FVA under 90% of max biomass production
    [minFlux(:, j), maxFlux(:, j)] = fluxVariabilityLLC(model, 90, 'max', model.rxns, 0, 1 - meth2test(j));
    fprintf('%s finished in %.2f sec.\n', methods{j}, toc(t0))
end

for i = 1:(numel(meth2test) - 1)
    for j = (i + 1):meth2test
        d = abs(minFlux(:, i) - minFlux(:, j));
        % check that all FVA values have same sign or the absolute differences are small
        assert(all(sign(minFlux(:, i)) == sign(minFlux(:, j)) | d < 1e-7))
        % check that for non-zero values the % differences are small
        f = minFlux(:, i) ~= 0 & minFlux(:, j) ~= 0;
        assert(max(d(f > 0) ./ f(f > 0)) < 1e-4)
    end
end
