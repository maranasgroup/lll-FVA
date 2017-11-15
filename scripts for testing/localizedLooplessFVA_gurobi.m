function [model,fvaRange, resInfo] = localizedLooplessFVA_gurobi(model, obj_vector, alpha, solverName, rxns, rxnLink, Nint, rxnInLoops, infoOnly, displayOn)
%%% This is only for comparing the localized loopless FVA and Fast-SNP
%%% preprocessed FVA. Please use the COBRA implementation
%%% fluxVariabilityLLC.m for general purpose.
% This function was modified from fastLooplessFVA by Saa and Nielsen 2016 to 
% have the same program structure for fair comparison. Only Gurobi can be used.
%
% Input:  
%     model:                 COBRA model structure
%
% Optional inputs:
%     obj_vector:            Objective vector to be maximized (default all zero)
%     alpha:                 % to optimality (default 0)
%     method:                'fast' or otherwise the traditional ll-FVA
%     solverName:            solver, only gurobi is implemented
%     rxns:                  reactions for FVA
%     rxnLink:               true to find reactions connected by EFMs in
%                            loops directly. Or provide the n-by-n matrix
%                            for the connection between reactions derived from EFMs directly
%     Nint:                  Provide null-space matrix for reactions in
%                            loops directly. Must be provided together with rxnInLoops
%     rxnInLoops:            n-by-2 matrix. rxnInLoops(j, 1) = true if the
%                            reverse direction of rxn j is in loop. 
%                            rxnInLoops(j, 2) = true if the forward
%                            direction of rxn j is in loop.
%     InfoOnly:              true to retrieve information on loop
%                            structures only, not performing FVA
%     displayOn:             true to print intermediate steps
%
% Outputs: 
%     model:                 structure with new bounds and directionalities, ll-FVA ranges
%     fvaRange:              FVA results
%     resInfo:               structure containing time cost and loop structure info.

% Check inputs
if ~exist('displayOn', 'var')
	displayOn = 0;
end
milpOutput = 0;

if nargin<1
    disp('Not enough input arguments.'); return;
elseif nargin < 2 || isempty(obj_vector)
    obj_vector = zeros(size(model.S, 2), 1);           % No obj fxn defined
    alpha  = 0;                % %-optimality
    if exist('gurobi','file')  % Check available solvers in this order
        solverName = 'gurobi';
    elseif exist('cplexlp','file')
        solverName = 'cplex';
    else
        addpath('glpkmex');    % We add glpk to the path
        solverName = 'glpk';
    end
elseif nargin < 3 || isempty(alpha)
    alpha  = 0;                % %-optimality
    if exist('gurobi','file')  % Check available solvers in this order
        solverName = 'gurobi';
    elseif exist('cplexlp','file')
        solverName = 'cplex';
    else
        addpath('glpkmex');    % We add glpk to the path
        solverName = 'glpk';
    end
elseif nargin < 4 || isempty(solverName)
    if exist('gurobi','file')  % Check available solvers in this order
        solverName = 'gurobi';
    elseif exist('cplexlp','file')
        solverName = 'cplex';
    else
        addpath('glpkmex');    % We add glpk to the path
        solverName = 'glpk';
    end
end
solverName = lower(solverName);
if nargin < 5 || isempty(rxns)
    rxns = 1:numel(model.rxns);
elseif ischar(rxns) && strcmp(rxns, 'info')
    rxns = [];
elseif iscell(rxns) || ischar(rxns)
    rxns = findRxnIDs(model, rxns);
elseif ~isnumeric(rxns)
    error('Incorrect rxns input!')
end

% get info only without performing FVA
if ~exist('infoOnly', 'var') || isempty(infoOnly)
    infoOnly = false;
end
% find rxn linkage by finding EFMs or not
useRxnLink = false;
if exist('rxnLink', 'var') && (~isscalar(rxnLink) || rxnLink)
    p = which('CalculateFluxModes');
    if ~isempty(p)
        useRxnLink = true;
    else
        fprintf('EFMtool not in Matlab path. Unable to calculate EFMs. Use nullspace information only.\n')
    end
end

% Assign parameters to the appropiate solver
if strcmp(solverName,'gurobi')
    options.outputflag      = 0;
    options.OptimalityTol   = 1e-9;
    options.FeasibilityTol  = 1e-9;
    options.IntFeasTol      = 1e-9;
    options.MIPGapAbs       = 1e-9;
    options.MIPGap          = 1e-9;
    options.Method = -1;
%     options.TimeLimit       = 5*60;    % 10 min time limit (uncomment if desired)
elseif strcmp(solverName,'cplex')
    options.Display         = 'off';
    options.TolFun          = 1e-9;
    options.TolRLPFun       = 1e-9; %1e-6;
    options.TolXInteger     = 1e-9; %1e-4;
    %     options.MaxTime         = 10*60;
else
    options.tolint          = 1e-9; %1e-4;
    options.tolobj          = 1e-9; %1e-6;
    options.tolbnd          = 1e-9; %1e-6;
    options.toldj           = 1e-9; %1e-6;
    %     options.tmlim           = 10*60;
end
tol = 1e-8;


%% Find minimal null-space and determine loop structures
% timer
t0 = cputime;
tNS = tic;
disp('Find minimal null-space...');
if ~exist('Nint', 'var') || ~exist('rxnInLoops', 'var') || isempty(Nint) || isempty(rxnInLoops)
    [rxnInLoops, Nint, loopInfo] = findMinNull_gurobi(model, 1, options);
    ex = sum(model.S ~= 0, 1) <= 1;
    rxnInLoops = rxnInLoops([find(~ex), find(ex)], :);
    Nint = Nint(~ex, :);
end
[model,rxnOrder] = parseInternalRxns_mod(model);
% record the time
tNS = toc(tNS);
cpuNS = cputime - t0;

[~, model.rxnOrderRev] = sort(rxnOrder);
rxns = model.rxnOrderRev(rxns);
solverName = lower(solverName);

%Initialize FVA ranges
[lpTime, milpTime, lpCPU, milpCPU, fvaRange, fvaObjBound] = deal(zeros(numel(rxns),2));
[mipStart, optiSol, presolveOff] = deal(false(numel(rxns),2));
nBinVar = zeros(numel(rxns), 1);

% check if the FVA must go completely loopless
t = tic;
tCPU = cputime;
obj_vector = obj_vector(:);

% find reactions in loops that are connected in the minimal null-space
conComp = -ones(size(model.S, 2), 1);
conComp(sum(model.S ~= 0, 1) > 1) = connectVars(Nint);

% determine the set of reactions for which LLCs are always required
% condition I in the paper (assume maximization of obj_vector)
cond1 = rxnInLoops(:, 2) & obj_vector > 0;
% condition II in the paper
cond2 = rxnInLoops(:, 1) & obj_vector < 0;
% condition III in the paper
cond3 = (model.lb > 0 & rxnInLoops(:, 2)) | (model.ub < 0 & rxnInLoops(:, 1));
rxnInLoopsAlwaysOn = cond1 | cond2 | cond3;
% LLCs are always required if the set is non-empty
alwaysLLC = any(rxnInLoopsAlwaysOn);
% the corresponding set of reactions in the same connected components as
% the always-on reactions
conCompAlwaysOn = false(max(conComp), 1);
conCompAlwaysOn(conComp(rxnInLoopsAlwaysOn)) = true;
loopPreprocessTime = toc(t);
loopPreprocessCPU = cputime - tCPU;
fprintf('Reactions in internal nullspace can be divided into %d connected components.\n', max(conComp))

if ~alwaysLLC
    % Get an LP model if LLCs are not always required, i.e., some problems
    % can be solved as LPs
    modelAllowLoops = allowLoopsLP(model, lower(solverName));
    objSolver = struct('gurobi', 'obj', 'cplex', 'f', 'glpk', 'c');
    modelAllowLoops.(objSolver.(lower(solverName)))(:) = 0;
end

% Find reactions connected by EFMs directly
if useRxnLink
    loadRxnLink = false;
    if ischar(rxnLink) && exist(rxnLink, 'file')
        rxnLinkLoad = rxnLink;
        rxnLink = load(rxnLinkLoad, 'rxnLink');
        rxnLink = rxnLink.rxnLink;
        if size(rxnLink, 1) == size(model.S, 2) && size(rxnLink, 2) == size(model.S, 2)
            loadRxnLink = true;
            rxnLink = [];
        end
    end
    if loadRxnLink || (ismatrix(rxnLink) && ...
            size(rxnLink, 1) == size(model.S, 2) && size(rxnLink, 2) == size(model.S, 2))
        rxnLinkTime = 0;
        rxnLinkCPU = 0;
    else
        % find rxnLink
        t = tic;
        tCPU = cputime;
        rxnLink = getRxnLink(model, Nint, conComp);
        rxnLinkTime = toc(t);
        rxnLinkCPU = cputime - tCPU;
    end
    fprintf('Use rxnLink matrix\n')
end
options.outputflag      = milpOutput;
options.DisplayInterval = 5;
% get an initial feasible loopless solution
model2 = struct();
[m, n] = size(model.S);
model2.A = [model.S,   sparse(m, n); ...                 % Sv        = 0
            speye(n), -speye(n); ...     %  v - |v| <= 0
           -speye(n), -speye(n)];        % -v - |v| <= 0
model2.rhs = zeros(m + n * 2, 1);
model2.lb = [model.lb; zeros(n, 1)];
model2.ub = [model.ub; max(abs(model.lb), abs(model.ub))];
model2.obj = [zeros(n, 1); ones(n, 1)];
model2.sense = char(['=' * ones(1, m), '<' * ones(1, n * 2)]);
sFeas = gurobi(model2, options);
x = sFeas.x(1:n);
clear sFeas model2

%% Set-up problem (GUROBI)
if ~isempty(rxns) && ~infoOnly
    if strcmp(solverName,'gurobi')
        
        model = looplessStructureMILP_GUROBI_mod(model,Nint);
        [rxnInLoop4L, rxnInLoop4U] = deal(1:2);
        rhs0 = model.rhs;
        fobj  = zeros(size(model.obj));
        
        % (1) Run initial FBA problem if indicated
        if any(obj_vector)
            model.obj(obj_vector(rxnOrder)~=0) = abs(obj_vector(obj_vector~=0));
            if sum(obj_vector)>0
                model.modelsense = 'max';
            else
                model.modelsense = 'min';
            end
            
            sol = gurobi(model,options);
            
            % Fix the value to alpha% of the optimal
            if sum(obj_vector)>0
                model.lb(model.obj~=0) = sol.objval*alpha;
                model.ub(model.obj~=0) = sol.objval;
            else
                model.lb(model.obj~=0) = sol.objval*(2-alpha); % This is in the case the optimal val is negative
                model.ub(model.obj~=0) = sol.objval;
            end
        end
        
        % Get an LP model if completely loopless FVA is not required
        if ~alwaysLLC
            if any(obj_vector)
                % Fix the value to alpha% of the optimal
                if sum(obj_vector)>0
                    modelAllowLoops.lb(model.obj(1:model.numRxns) ~= 0) = sol.objval*alpha;
                    modelAllowLoops.ub(model.obj(1:model.numRxns) ~= 0) = sol.objval;
                else
                    modelAllowLoops.lb(model.obj(1:model.numRxns) ~= 0) = sol.objval*(2-alpha); % This is in the case the optimal val is negative
                    modelAllowLoops.ub(model.obj(1:model.numRxns) ~= 0) = sol.objval;
                end
            end
            modelAllowLoops.obj = fobj(1:model.numRxns);
        end
        
        model.obj = fobj;
        
        
        % (2) ll-FVA optimization loop
        for i0 = 1:numel(rxns) 
            %%
            i = rxns(i0);
            % Restart objective
            if any(~rxnInLoops(i, :))
                modelAllowLoops.obj(:) = 0;
            end
            % if the reaction for FVA is in loops
            if alwaysLLC || any(rxnInLoops(i, :))
                % reset objective, RHS, lb and ub
                model.obj(:) = 0;
                model.rhs = rhs0;
                model.ub(model.milpInfo.var.z) = 1;
                model.ub(model.milpInfo.var.g) = model.milpInfo.K;
                model.lb(model.milpInfo.var.g) = -model.milpInfo.K;
                % if using rxn connections from EFMs
                if useRxnLink
                    % if load given file for the connection matrix
                    if loadRxnLink
                        rxnLink = load(rxnLinkLoad, 'rxnLink');
                        rxnLink = rxnLink.rxnLink;
                    end
                    % relax constraints not connected to the reaction for FVA 
                    % except reactions required to be always constrained
                    bigM = inf;
                    id = rxnLink(i, :)' == 0 & any(rxnInLoops(:, rxnInLoop4U), 2) & ~rxnInLoopsAlwaysOn;
                    model.rhs(model.milpInfo.con.vU(model.milpInfo.rxnFwdLoopIds(id))) = bigM;
                    model.rhs(model.milpInfo.con.gU(model.milpInfo.rxnFwdLoopIds(id))) = bigM;
                    id = rxnLink(i, :)' == 0 & any(rxnInLoops(:, rxnInLoop4L), 2) & ~rxnInLoopsAlwaysOn;
                    model.rhs(model.milpInfo.con.vL(model.milpInfo.rxnRevLoopIds(id))) = bigM;
                    model.rhs(model.milpInfo.con.gL(model.milpInfo.rxnRevLoopIds(id))) = bigM;
                    % pre-determine variables not connected to the reaction for FVA
                    % except reactions required to be always constrained
                    model.lb(model.milpInfo.var.g(model.milpInfo.rxnLoopIds(conComp ~= conComp(i) & conComp > 0 & ~rxnInLoopsAlwaysOn))) = 0;
                    model.ub(model.milpInfo.var.g(model.milpInfo.rxnLoopIds(conComp ~= conComp(i) & conComp > 0 & ~rxnInLoopsAlwaysOn))) = 0;
                    model.ub(model.milpInfo.var.z(model.milpInfo.rxnLoopIds(rxnLink(i, :)' == 0 & any(rxnInLoops, 2) & ~rxnInLoopsAlwaysOn))) = 0;
                    nBinVar(i0) = sum(rxnLink(i, :)' | rxnInLoopsAlwaysOn);
                    if loadRxnLink
                        clear rxnLink
                    end
                else
                    conCompOnCur = conCompAlwaysOn;
                    conCompOnCur(conComp(i)) = true;
                    for jCon = 1:numel(conCompOnCur)
                        if ~conCompOnCur(jCon)
                            % relax the constraint
                            model.rhs(model.milpInfo.con.vU(model.milpInfo.rxnFwdLoopIds(conComp == jCon & any(rxnInLoops(:, rxnInLoop4U), 2)))) = inf;
                            model.rhs(model.milpInfo.con.gU(model.milpInfo.rxnFwdLoopIds(conComp == jCon & any(rxnInLoops(:, rxnInLoop4U), 2)))) = inf;
                            model.rhs(model.milpInfo.con.vL(model.milpInfo.rxnRevLoopIds(conComp == jCon & any(rxnInLoops(:, rxnInLoop4L), 2)))) = inf;
                            model.rhs(model.milpInfo.con.gL(model.milpInfo.rxnRevLoopIds(conComp == jCon & any(rxnInLoops(:, rxnInLoop4L), 2)))) = inf;
                            % fix variables not involved
                            model.lb(model.milpInfo.var.g(model.milpInfo.rxnLoopIds(conComp == jCon))) = 0;
                            model.ub(model.milpInfo.var.g(model.milpInfo.rxnLoopIds(conComp == jCon))) = 0;
                            model.ub(model.milpInfo.var.z(model.milpInfo.rxnLoopIds(conComp == jCon))) = 0;
                        end
                    end
                    nBinVar(i0) = sum(ismember(conComp, find(conCompOnCur)));
                end
            end
            if displayOn
				fprintf('Running %d / %d.\t%d binary vars.\t%04d-%02d-%02d %02d:%02d:%02.0f\n', i0, numel(rxns), nBinVar(i0), clock);
            end
            % Minimization problem
            if ~alwaysLLC && ~rxnInLoops(i, 1)
                % simple FVA
                modelAllowLoops.obj(i) = 1;
                modelAllowLoops.modelsense = 'min';
                t = tic;
                t1 = cputime;
                sol           = gurobi(modelAllowLoops,options);
                lpTime(i0, 1) = toc(t);
                lpCPU(i0, 1) = cputime - t1;
                fvaObjBound(i0,1) = sol.objval*(abs(sol.objval)>tol);
            else
                % loopless FVA if the objective vector contains rxns in loops
                % or the current rxn is in loops
                model.obj(i) = 1;
                model.modelsense = 'min';
                t = tic;
                t1 = cputime;
                sol           = gurobi(model,options);
                milpTime(i0, 1) = toc(t);
                milpCPU(i0, 1) = cputime - t1;
                if strcmpi(sol.status, 'INFEASIBLE')
                    % supply MipStart if infeasible
                    fprintf('MipStart needed.\n')
                    mipStart(i0, 1) = true;
                    xFeas = getMipStart_gurobi(x, Nint, conComp, model.ub(model.milpInfo.var.z), options);
                    model.start = xFeas;
                    t = tic;
                    t1 = cputime;
                    sol = gurobi(model, options);
                    milpTime(i0, 1) = toc(t);
                    milpCPU(i0, 1) = cputime - t1;
                    % turn off presolve if still infeasible
                    if strcmpi(sol.status, 'INFEASIBLE')
                        presolveOff(i0, 1) = true;
						options.Presolve = 0;
						t = tic;
                        t1 = cputime;
						sol = gurobi(model, options);
						milpTime(i0, 1) = toc(t);
                        milpCPU(i0, 1) = cputime - t1;
						options.Presolve = -1;
                    end
                    model = rmfield(model, 'start');
                end
                if strcmpi(sol.status, 'OPTIMAL')
					optiSol(i0, 1) = true;
                end
                fvaObjBound(i0, 1) = sol.objbound * (abs(sol.objbound) > tol);
            end
            fvaRange(i0,1) = sol.objval*(abs(sol.objval)>tol);
            
            if displayOn
				fprintf('    min: stat %s, obj %.4f, bound %.4f, %.2f sec\n', ...
					sol.status, fvaRange(i0, 1), fvaObjBound(i0, 1), toc(t))
            end
            pause(1e-3)
            % Maximization problem
            if ~alwaysLLC && ~rxnInLoops(i, 2)
                % simple FVA
                modelAllowLoops.obj(i) = 1;
                modelAllowLoops.modelsense = 'max';
                t = tic;
                t1 = cputime;
                sol           = gurobi(modelAllowLoops,options);
                lpTime(i0, 2) = toc(t);
                lpCPU(i0, 2) = cputime - t1;
                fvaObjBound(i0,2) = sol.objval*(abs(sol.objval)>tol);
            else
                model.obj(i) = 1;
                model.modelsense = 'max';
                t = tic;
                t1 = cputime;
                sol           = gurobi(model,options);
                milpTime(i0, 2) = toc(t);
                milpCPU(i0, 2) = cputime - t1;
                % supply MipStart if infeasible
                if strcmpi(sol.status, 'INFEASIBLE')
                    fprintf('MipStart needed.\n')
                    mipStart(i0, 2) = true;
                    xFeas = getMipStart_gurobi(x, Nint, conComp, model.ub(model.milpInfo.var.z), options);
                    model.start = xFeas;
                    t = tic;
                    t1 = cputime;
                    sol = gurobi(model, options);
                    milpTime(i0, 2) = toc(t);
                    milpCPU(i0, 2) = cputime - t1;
                    % turn off presolve if still infeasible
                    if strcmpi(sol.status, 'INFEASIBLE')
                        presolveOff(i0, 2) = true;
						options.Presolve = 0;
						t = tic;
                        t1 = cputime;
						sol = gurobi(model, options);
						milpTime(i0, 2) = toc(t);
                        milpCPU(i0, 2) = cputime - t1;
						options.Presolve = -1;
                    end
                    model = rmfield(model, 'start');
                end
                if strcmpi(sol.status, 'OPTIMAL')
					optiSol(i0, 2) = true;
                end
                fvaObjBound(i0, 2) = sol.objbound * (abs(sol.objbound) > tol);
            end
            %%
            fvaRange(i0,2) = sol.objval*(abs(sol.objval)>tol);
            
            if displayOn
				fprintf('    max: stat %s, obj %.4f, bound %.4f, %.2f sec\n', ...
					sol.status, fvaRange(i0, 2), fvaObjBound(i0, 2), toc(t))
            end
            pause(1e-3)
            
        end
        
        %% Set-up problem (CPLEX)
    elseif strcmp(solverName,'cplex')
        model = looplessStructureMILP_CPLEX(model,Nint);
        model.options = options;
        
        fobj  = zeros(size(model.f));
        
        % Run initial FBA problem if indicated
        if any(obj_vector)
            model.f(obj_vector(rxnOrder)~=0) = -obj_vector(obj_vector~=0); % deals with either max or min
            [~,fval] = cplexmilp(model);
            
            % Fix the value to alpha% of the optimal
            if sum(obj_vector)>0
                model.lb(model.f~=0) = -fval*alpha;
                model.ub(model.f~=0) = -fval;
            else
                model.lb(model.f~=0) = -fval*(2-alpha);
                model.ub(model.f~=0) = -fval;
            end
        end
        
        % Get an LP model if completely loopless FVA is not required
        if ~alwaysLLC
            if any(obj_vector)
                % Fix the value to alpha% of the optimal
                if sum(obj_vector)>0
                    modelAllowLoops.lb(model.f(1:model.numRxns) ~= 0) = fval*alpha;
                    modelAllowLoops.ub(model.f(1:model.numRxns) ~= 0) = fval;
                else
                    modelAllowLoops.lb(model.f(1:model.numRxns) ~= 0) = fval*(2-alpha); % This is in the case the optimal val is negative
                    modelAllowLoops.ub(model.f(1:model.numRxns) ~= 0) = fval;
                end
            end
            modelAllowLoops.f = fobj(1:model.numRxns);
            modelAllowLoops.options = options;
        end
        
        model.f = fobj;
        
        % (2) ll-FVA optimization loop
        for i = 1:model.numRxns
            
            % Minimization problem
            if ~alwaysLLC && ~rxnInLoops(i, 1)
                % simple FVA
                modelAllowLoops.f(i) = 1;
                [~,fval]      = cplexlp(modelAllowLoops);
            else
                % loopless FVA
                model.f(i)    = 1;
                [~,fval]      = cplexmilp(model);
            end
            fvaRange(i,1) = fval*(abs(fval)>tol);
            
            % Maximization problem
            if ~alwaysLLC && ~rxnInLoops(i, 2)
                modelAllowLoops.f(i) = -1;
                [~,fval]      = cplexlp(modelAllowLoops);
                % Restart objective
                modelAllowLoops.f = fobj(1:model.numRxns);
            else
                % loopless FVA
                model.f(i)    = -1;
                [~,fval]      = cplexmilp(model);
                % Restart objective
                model.f = fobj;
            end
            fvaRange(i,2) = -fval*(abs(fval)>tol);
            
        end
        
        %% Set-up problem (GLPK)
    elseif strcmp(solverName,'glpk')
        model = looplessStructureMILP_GLPK(model,Nint);
        
        fobj  = zeros(size(model.c));
        
        % Run initial FBA problem if indicated
        if any(obj_vector)
            model.c(obj_vector(rxnOrder)~=0) = abs(obj_vector(obj_vector~=0)); % deals with either max or min
            [~,fval] = glpk(model.c,model.A,model.b,model.lb,model.ub,model.ctype,model.vartype,model.sense,options);
            
            % Fix the value to alpha% of the optimal
            if sum(obj_vector)>0
                model.lb(model.c~=0) = fval*alpha;
                model.ub(model.c~=0) = fval;
            else
                model.lb(model.c~=0) = fval*(2-alpha);
                model.ub(model.c~=0) = fval;
            end
        end
        
        if ~alwaysLLC
            if any(obj_vector)
                % Fix the value to alpha% of the optimal
                if sum(obj_vector)>0
                    modelAllowLoops.lb(model.c(1:model.numRxns) ~= 0) = fval*alpha;
                    modelAllowLoops.ub(model.c(1:model.numRxns) ~= 0) = fval;
                else
                    modelAllowLoops.lb(model.c(1:model.numRxns) ~= 0) = fval*(2-alpha); % This is in the case the optimal val is negative
                    modelAllowLoops.ub(model.c(1:model.numRxns) ~= 0) = fval;
                end
            end
            modelAllowLoops.c = fobj(1:model.numRxns);
        end
        
        model.c = fobj;
        
        % (2) ll-FVA optimization loop
        for i = 1:model.numRxns
            model.c(i)    = 1;
            
            % Minimization problem
            % Minimization problem
            if ~alwaysLLC && ~rxnInLoops(i, 1)
                % simple FVA
                modelAllowLoops.sense = 1;
                [~,fval]      = glpk(modelAllowLoops.c, modelAllowLoops.A, ...
                    modelAllowLoops.b, modelAllowLoops.lb, modelAllowLoops.ub, ...
                    modelAllowLoops.ctype, modelAllowLoops.vartype, ...
                    modelAllowLoops.sense, options);
            else
                model.sense   = 1;
                [~,fval]      = glpk(model.c,model.A,model.b,model.lb,model.ub,model.ctype,model.vartype,model.sense,options);
            end
            
            fvaRange(i,1) = fval*(abs(fval)>tol);
            
            % Maximization problem
            if ~alwaysLLC && ~rxnInLoops(i, 2)
                % simple FVA
                modelAllowLoops.sense = -1;
                [~,fval]      = glpk(modelAllowLoops.c, modelAllowLoops.A, ...
                    modelAllowLoops.b, modelAllowLoops.lb, modelAllowLoops.ub, ...
                    modelAllowLoops.ctype, modelAllowLoops.vartype, ...
                    modelAllowLoops.sense, options);
                % Restart objective
                modelAllowLoops.c = fobj(1:model.numRxns);
            else
                model.sense   = -1;
                [~,fval]      = glpk(model.c,model.A,model.b,model.lb,model.ub,model.ctype,model.vartype,model.sense,options);
                % Restart objective
                model.c = fobj;
            end
            
            fvaRange(i,2) = fval*(abs(fval)>tol);
        end
    end
    
    % Save ll-FVA time
    elapsedTime = (cputime-t0);
    
    % Find zero reactions
    zeroRxns = [];%rxns(~any(fvaRange, 2));  % abs(fvaRange(:,1))+abs(fvaRange(:,2))==0);
    if isempty(zeroRxns)
        
        % Reset objective function
        model.obj = fobj;
        
        % Redefine model boundaries
        model.lb(rxns) = fvaRange(:,1);
        model.ub(rxns) = fvaRange(:,2);
    else
        % Remove zero rxns
        modelTemp.S = model.S;
        modelTemp.S(:,zeroRxns) = [];
        modelTemp.c = zeros(size(modelTemp.S,2),1);
        Nint        = Nint(~ismember(model.internal,zeroRxns),:);
        
        % Find orphan metabolites
        orphanMets = find(sum(full(abs(modelTemp.S)),2)==0);
        modelTemp.S(orphanMets,:) = [];
        
        % Re-define other quantities
        modelTemp.b                    = zeros(size(modelTemp.S,1),1);
        [modelTemp.lb, modelTemp.ub]   = deal(model.lb(1:model.numRxns), model.ub(1:model.numRxns));
        modelTemp.lb(rxns)             = fvaRange(:,1);
        modelTemp.lb(zeroRxns)         = [];
        modelTemp.ub(rxns)             = fvaRange(:,2);
        modelTemp.ub(zeroRxns)         = [];
        modelTemp.rxns                 = model.rxns;
        modelTemp.rxns(zeroRxns)       = [];
        modelTemp.rxnNames             = model.rxnNames;
        modelTemp.rxnNames(zeroRxns)   = [];
        modelTemp.mets                 = model.mets;
        modelTemp.mets(orphanMets)     = [];
        modelTemp.metNames             = model.metNames;
        modelTemp.metNames(orphanMets) = [];
        modelTemp.description          = model.description;
        
        % Re-format the problem
        if strcmp(solverName,'gurobi')
			model = looplessStructureMILP_GUROBI_mod(parseInternalRxns_mod(modelTemp),Nint);
        elseif strcmp(solverName,'cplex')
            model = looplessStructureMILP_CPLEX(parseInternalRxns_mod(modelTemp),Nint);
        elseif strcmp(solverName,'glpk')
            model = looplessStructureMILP_GLPK(parseInternalRxns_mod(modelTemp),Nint);
        end
    end
    
    % Assign reversibilities to thermodynamically allowable rxns
    model.rev         = (model.lb(1:model.numRxns)<0).*(model.ub(1:model.numRxns)>0);
    disp(['Done. Reactions removed ',num2str(length(zeroRxns))]);
else
    elapsedTime = cputime - t0;
    fprintf('Get information on loops only.\n')
end

model.prepTime    = elapsedTime;
model.rxnOrder = rxnOrder;

resInfo = struct();
if isempty(rxns) || infoOnly
    resInfo.Nint = Nint;
end
resInfo.rxns = rxns;
resInfo.rxns2 = model.rxns(rxns);
resInfo.rxnInLoops = rxnInLoops;
resInfo.conComp = conComp;
resInfo.nBinVar = nBinVar;
resInfo.nsTime = tNS - loopInfo.loopPreprocessTime;
resInfo.nsCPU = cpuNS - loopInfo.loopPreprocessCPU;
resInfo.loopPreprocessTime = loopPreprocessTime + loopInfo.loopPreprocessTime;
resInfo.loopPreprocessCPU = loopPreprocessCPU + loopInfo.loopPreprocessCPU;
resInfo.lpTime = lpTime;
resInfo.milpTime = milpTime;
resInfo.lpCPU = lpCPU;
resInfo.milpCPU = milpCPU;
resInfo.mipStart = mipStart;
resInfo.optiSol = optiSol;
resInfo.presolveOff = presolveOff;
resInfo.loopInfo = loopInfo;
if useRxnLink
    resInfo.rxnLink = rxnLink;
    resInfo.rxnLinkTime = rxnLinkTime;
    resInfo.rxnLinkCPU = rxnLinkCPU;
else
    resInfo.rxnLinkTime = 0;
    resInfo.rxnLinkCPU = 0;
end
resInfo.fvaObjBound = fvaObjBound;


function modelAllowLoops = allowLoopsLP(model, solver)
modelAllowLoops = struct();
[nMets, nRxns] = size(model.S);
switch solver
    case 'gurobi'
        modelAllowLoops.A = sparse(model.S);
        modelAllowLoops.rhs = zeros(nMets, 1);
        modelAllowLoops.sense = char('=' * ones(1, nMets));
        modelAllowLoops.vtype = char('C' * ones(1, nRxns));
        modelAllowLoops.obj = zeros(nRxns, 1);
    case 'cplex'
        modelAllowLoops.Aeq = model.S;
        modelAllowLoops.beq = zeros(nMets, 1);
        modelAllowLoops.Aineq = [];
        modelAllowLoops.bineq = [];
        modelAllowLoops.f = zeros(nRxns, 1);
    case 'glpk'
        modelAllowLoops.A = model.S;
        modelAllowLoops.b = zeros(nMets, 1);
        modelAllowLoops.ctype = char('S' * ones(1, nMets));
        modelAllowLoops.vartype = char('C' * ones(1, nRxns));
        modelAllowLoops.c = zeros(nRxns, 1);
end
modelAllowLoops.lb = model.lb(1:nRxns);
modelAllowLoops.ub = model.ub(1:nRxns);
