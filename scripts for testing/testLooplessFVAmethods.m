function testLooplessFVAmethods(testModel, testVector, testMethod, saveName, modelPath)
% Test the 7 models tested in the paper of Fast-SNP (Saa and Nielson, 2016)

% set the scope of the test here
% set the data to plot here
if nargin < 1 || isempty(testModel)
    % models to be tested (1 to 7)
    % (The first three models can be finised within 2~3 mins for each run.
    %  The last four models can take 10~60 mins for ll-FVA. But still fast for
    %  LLC methods)
    testModel = 1:3;
end
if nargin < 2 || isempty(testVector)
    % number of tests to perform
    testVector = 1:3;
end
if nargin < 3 || isempty(testMethod)
    % methods to be tested (1: 'll-FVA', 2: 'fastSNP', 3: 'LLC-NS', 4: 'LLC-EFM');
    testMethod = 1:4;
end
if nargin < 4 || isempty(saveName)
    % path for saving data files
    saveName = [pwd filesep 'test_results' filesep 'testResult_'];
end
if nargin < 5 || isempty(modelPath)
    modelPath = [pwd filesep 'test_models'];
end

allMethod = {'ll-FVA', 'fastSNP', 'NS-LLC', 'EFM-LLC'};
method2Test = allMethod(testMethod);  % methods to be tested

%% get the model list
nTest = numel(testVector);
d = dir(modelPath);
modelList = {};
modelSize = [];
for j = 1:numel(d)
    if ~strncmp(d(j).name(1), '.', 1) && ~strncmp(d(j).name, 'model', 5)
        modelList(end + 1) = {[modelPath, filesep, d(j).name]};
        model = load(modelList{end});
        modelSize = [modelSize; size(model.model.S)];
    end
end
% start from small models to large models
[~, ind] = sort(modelSize(:, 1) .* modelSize(:, 2));
modelList = modelList(ind)';
modelSize = modelSize(ind, :);
modelName = regexp(modelList, ['[^\' filesep ']+$'], 'match', 'once');

disp([{'model to test', '#mets', '#rxns'}; modelName num2cell(modelSize)])

solver = {'gurobi'};

saveDir = strsplit(saveName, filesep);
saveDir = strjoin(saveDir(1:end - 1), filesep);
if ~isempty(saveDir) && ~exist(saveDir, 'dir')
    mkdir(saveDir)
end

%% loop for each model
for j = testModel(:)'
    
    model = load(modelList{j});
    fprintf('\nTest model %s\n', modelName{j});
    model = model.model;
    
    % Run without any objective. Just get the bounds for the models
    model.c(:) = 0;
    rxns = model.rxns;
    for meth = 1:numel(method2Test)
        for jTest = testVector(:)'
            fprintf('Test %d Method %s\n', jTest, method2Test{meth})
            saveNameJ = [saveName 't' num2str(jTest)];
            if ~exist([saveNameJ '_m' num2str(j) '_' method2Test{meth} '.mat'], 'file')
                clear modelRes fvaRange resInfo
                t0 = tic;
                t1 = cputime;
                switch method2Test{meth}
                    case 'll-FVA'  % original loopless FVA
                        [modelRes, fvaRange, resInfo] = fastLooplessFVA_mod(model, [], 0, 'original', 'gurobi', rxns);
                    case 'fastSNP'  % loopless FVA with Fast-SNP preprocessing
                        [modelRes, fvaRange, resInfo] = fastLooplessFVA_mod(model, [], 0, 'fast', 'gurobi', rxns);
                    case 'NS-LLC'  % null-space-based localized loopless FVA
                        [modelRes, fvaRange, resInfo] = localizedLooplessFVA_gurobi(model, [], 0, 'gurobi', rxns, false);
                    case 'EFM-LLC'  % EFM-based localized loopless FVA
                        [modelRes, fvaRange, resInfo] = localizedLooplessFVA_gurobi(model, [], 0, 'gurobi', rxns, true);
                    otherwise
                        error('Method unknown!')
                end
                wallTime = toc(t0);
                cpuTime = cputime - t1;
                pause(1e-3)
                nRxnCyc = sum(any(modelRes.Nint, 1));
                fprintf('%s finished in %.4f sec.\n', method2Test{meth}, wallTime)
                fprintf('%04d-%02d-%02d %02d:%02d:%02.0f\n\n', clock)
                save([saveNameJ '_m' num2str(j) '_' method2Test{meth} '.mat'], 'fvaRange', ...
                    'cpuTime', 'nRxnCyc', 'modelList', 'modelName', ...
                    'method2Test','meth', 'wallTime', 'resInfo', 'rxns', '-v7.3');
                % ensure that the same set of reactions was checked
                assert(isempty(setdiff(modelRes.rxns(resInfo.rxns),rxns)) & numel(modelRes.rxns(resInfo.rxns)) == numel(rxns))
            end
        end
    end
end
%%
if 1
%%
    clear data
    fprintf('\nMax diff between methods:\n')
    fprintf('#test  model             ');
    [maxDiff, maxDiffP] = deal(zeros(nTest, numel(method2Test), numel(testModel)));

    for jMethod = 1:numel(method2Test)
        fprintf('%-10s  ', method2Test{jMethod})
    end
    fprintf('\n')
    for jTest = testVector(:)'
        saveNameJ = [saveName 't' num2str(jTest)];
        for j = testModel(:)' 
            fprintf('%5d  %-15s  ', jTest, modelName{j})
            for jMethod = 1:numel(method2Test)
                data(jMethod) = load([saveNameJ '_m' num2str(j) '_' method2Test{jMethod} '.mat']);
                maxDiff(jTest, jMethod, j) = max(max(abs(data(jMethod).fvaRange - data(1).fvaRange)));
                maxV = max(abs(data(jMethod).fvaRange), abs(data(1).fvaRange));
                maxDiffP(jTest, jMethod, j) = max(max( abs(data(jMethod).fvaRange(maxV ~= 0) - data(1).fvaRange(maxV ~= 0)) ./ maxV(maxV ~= 0)));
                fprintf('%-10.4e  ', maxDiff(jTest, jMethod, j))
            end
            fprintf('\n')
        end
    end
    fprintf('Max percentage diff:\n')
    for jTest = testVector(:)'
        for j = testModel(:)' 
            fprintf('%5d  %-15s  ', jTest, modelName{j})
            for jMethod = 1:numel(method2Test)
                fprintf('%-10.4e  ', maxDiffP(jTest, jMethod, j))
            end
            fprintf('\n')
        end
    end

end
