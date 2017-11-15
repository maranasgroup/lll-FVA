function [model,fvaRange, resInfo] = fastLooplessFVA_mod(model, obj_vector, alpha, method, solverName, rxns, Nint, displayOn)
% Performs conventional or fast ll-FVA. In the latter case, Fast-SNP is
% employed to find a suitable basis for Nint
% Inputs:  model structure, obj_vector (optional), alpha (% to optimality, optional)
% Outputs: model structure with new bounds and directionalities, ll-FVA
%          ranges
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2016 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Improved parsing of matrices and take additional inputs. Only applicable to using Gurobi
% Added inputs:
%    rxns: reactions for which loopless FVA will be performed
%    Nint: nulls-pace matrix for internal reactions if already known
%    displayOn: true to print intermediate steps
%
% Added output:
%    resInfo: structure containing time cost and other info.
% (commented starting with '%^')
%%%%% Siu Hung Joshua Chan 2017

% Check inputs

%^ print intermediate steps
if ~exist('displayOn', 'var')
	displayOn = 0;
end

if nargin<1
    disp('Not enough input arguments.'); return;
end
if nargin < 2 || isempty(obj_vector)
    obj_vector = zeros(size(model.S, 2), 1);  % No obj fxn defined
end
if nargin < 3 || isempty(alpha)
    alpha  = 0;                % %-optimality
end
if nargin < 4  || isempty(method)
    method = 'fast';
end
if nargin < 5 || isempty(solverName)
    if exist('gurobi','file')  % Check available solvers in this order
        solverName = 'gurobi';
    elseif exist('cplexlp','file')
        solverName = 'cplex';
    else
        addpath('glpkmex');    % We add glpk to the path
        solverName = 'glpk';
    end
end
%^ specify the set of reactions for FVA
if nargin < 6 || isempty(rxns)
    rxns = 1:numel(model.rxns);
elseif ischar(rxns) && strcmp(rxns, 'info')
    rxns = [];
elseif iscell(rxns) || ischar(rxns)
    rxns = findRxnIDs(model, rxns);
elseif ~isnumeric(rxns)
    error('Incorrect rxns input!')
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
    % options.Presolve = 0;
    % options.TimeLimit       = 10*60;    % 10 min time limit (uncomment if desired)
elseif strcmp(solverName,'cplex')
    options.Display         = 'off';
    options.TolFun          = 1e-9;
    options.TolRLPFun       = 1e-9;
    options.TolXInteger     = 1e-9;
    %     options.MaxTime         = 10*60;
else
    options.tolint          = 1e-4;
    options.tolobj          = 1e-6;
    options.tolbnd          = 1e-6;
    options.toldj           = 1e-6;
    %     options.tmlim           = 10*60;
end
tol = 1e-8;

% Main script
t0 = cputime;
tSNP = tic;

% Determine sparse basis for Nint
if strcmp('fast',method)
    disp('Performing ll-FVA using Fast-SNP...');
    if ~exist('Nint', 'var') || isempty(Nint)  %^ added to allow user-supplied Nint
		Nint  = fastSNP_mod(model,solverName);
	end
    [model,rxnOrder] = parseInternalRxns_mod(model);
else
    disp('Performing traditional ll-FVA...');
    [model,rxnOrder] = parseInternalRxns_mod(model);
    if ~exist('Nint', 'var') || isempty(Nint)  %^ added to allow user-supplied Nint
		Nint = getNullSpace(model.S(:, model.internal), 0);
	end
    Nint(abs(Nint)<tol) = 0; % Remove elements below numerical tolerance
end

%^ add timer
tSNP = toc(tSNP);
cpuSNP = cputime - t0;
%^ sort the order
[~, model.rxnOrderRev] = sort(rxnOrder);
rxns = model.rxnOrderRev(rxns);

% Initialize FVA ranges  %^ Added other variables
[lpTime, milpTime, lpCPU, milpCPU, fvaRange] = deal(zeros(numel(rxns),2));
[mipStart, optiSol, presolveOff] = deal(false(numel(rxns),2));
nBinVar = zeros(numel(rxns), 1);
fvaObjBound = zeros(size(fvaRange));

%^ rxnInLoops(j, 1) = true means the reverse direction of rxn j is in loops
%^ rxnInLoops(j, 2) = true means the forward direction of rxn j is in loops
rxnInLoops = [any(Nint, 2); false(size(model.S, 2) - size(Nint, 1), 1)];
nRxnInLoops = sum(any(rxnInLoops, 2));

%^ get an initial feasible loopless solution in case sometimes solvers with
%^ presolve turned on would result an infeasible problem
model2 = model;
model2.lb = model2.lb(1:size(model2.S, 2));
model2.ub = model2.ub(1:size(model2.S, 2));
model2.c = zeros(size(model2.S, 2), 1);
model2.b = zeros(size(model2.S, 1), 1);
sFeas = optimizeCbModel(model2, 'max', 'one');
xFeas = getMipStart_gurobi(sFeas.x, Nint, [], [], options);
clear sFeas model2


%% Set-up problem (GUROBI)
%^ only gurobi was edited. Other solvers were unchanged.
if strcmp(solverName,'gurobi')
    %^ use sparse matrix to increase parsing speed
    model = looplessStructureMILP_GUROBI_mod(model,Nint);
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
    model.obj = fobj;
    
    % (2) ll-FVA optimization loop
    for i0 = 1:numel(rxns)
		i = rxns(i0);
        model.obj(i) = 1;
        %^ record number of binary varibles in solving the MILPs
        nBinVar(i0) = nRxnInLoops;
        if displayOn  % print for each solve
            fprintf('Running %d / %d.\t%d binary vars.\t%04d-%02d-%02d %02d:%02d:%02.0f\n', i0, numel(rxns), nBinVar(i0), clock);
        end
        % Minimization problem
        model.modelsense = 'min';
        %^ timer
        t = tic;
        t1 = cputime;
        sol           = gurobi(model,options);
        %^ record the time
        milpTime(i0, 1) = toc(t);
        milpCPU(i0, 1) = cputime - t1;
        %^ handle infeasible solution due to presolve
        if strcmpi(sol.status, 'INFEASIBLE')
            fprintf('MipStart needed.\n')
            %^ give an initial solution
            mipStart(i0, 1) = true;
            model.start = xFeas;
            t = tic;
            t1 = cputime;
            sol = gurobi(model, options);
            milpTime(i0, 1) = toc(t);
            milpCPU(i0, 1) = cputime - t1;
            %^ if still not feasible, abort presolve
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
        %^ check if it is the optimal solution
        if strcmpi(sol.status, 'OPTIMAL')
			optiSol(i0, 1) = true;
        end
        
        fvaRange(i0,1) = sol.objval*(abs(sol.objval)>tol);
        fvaObjBound(i0, 1) = sol.objbound * (abs(sol.objbound) > tol);
        
        %^ solver's behavior was seen weird occasionally in the tested computer if not giving a short pause
        pause(1e-5)
        if displayOn
            fprintf('    min: stat %s, obj %.4f, bound %.4f, %.2f sec\n', ...
                sol.status, fvaRange(i0, 1), fvaObjBound(i0, 1), toc(t))
        end
        
        % Maximization problem
        model.modelsense = 'max';
        t = tic;
        t1 = cputime;
        sol           = gurobi(model,options);
        %^ record the time
        milpTime(i0, 2) = toc(t);
        milpCPU(i0, 2) = cputime - t1;
        %^ handle infeasible solution due to presolve
        if strcmpi(sol.status, 'INFEASIBLE')
            fprintf('MipStart needed.\n')
            %^ give an initial solution
            mipStart(i0, 2) = true;
            model.start = xFeas;
            t = tic;
            t1 = cputime;
            sol = gurobi(model, options);
            milpTime(i0, 2) = toc(t);
            milpCPU(i0, 2) = cputime - t1;
            %^ if still not feasible, abort presolve
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
        %^ check if it is the optimal solution
        if strcmpi(sol.status, 'OPTIMAL')
			optiSol(i0, 2) = true;
        end
        
        fvaRange(i0,2) = sol.objval*(abs(sol.objval)>tol);
        fvaObjBound(i0, 2) = sol.objbound * (abs(sol.objbound) > tol);
        % Restart objective
        model.obj = fobj;
        
        %^ solver's behavior was seen weird occasionally in the tested computer if not giving a short pause
        pause(1e-5)
        if displayOn
            fprintf('    max: stat %s, obj %.4f, bound %.4f, %.2f sec\n', ...
                sol.status, fvaRange(i0, 2), fvaObjBound(i0, 2), toc(t))
        end
    end
    
    %% Set-up problem (CPLEX)
elseif strcmp(solverName,'cplex')
    model = looplessStructureMILP_CPLEX(model,Nint);
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
    model.f = fobj;
    
    % (2) ll-FVA optimization loop
    for i = 1:model.numRxns
        
        % Minimization problem
        model.f(i)    = 1;
        [~,fval]      = cplexmilp(model);
        fvaRange(i,1) = fval*(abs(fval)>tol);
        
        % Maximization problem
        model.f(i)    = -1;
        [~,fval]      = cplexmilp(model);
        fvaRange(i,2) = -fval*(abs(fval)>tol);
        
        % Restart objective
        model.f = fobj;
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
    model.c = fobj;
    
    % (2) ll-FVA optimization loop
    for i = 1:model.numRxns
        model.c(i)    = 1;
        
        % Minimization problem
        model.sense   = 1;
        [~,fval]      = glpk(model.c,model.A,model.b,model.lb,model.ub,model.ctype,model.vartype,model.sense,options);
        fvaRange(i,1) = fval*(abs(fval)>tol);
        
        % Maximization problem
        model.sense   = -1;
        [~,fval]      = glpk(model.c,model.A,model.b,model.lb,model.ub,model.ctype,model.vartype,model.sense,options);
        fvaRange(i,2) = fval*(abs(fval)>tol);
        
        % Restart objective
        model.c = fobj;
    end
end

% Save ll-FVA time
elapsedTime = (cputime-t0);

% Find zero reactions
zeroRxns = [];%find(abs(fvaRange(:,1))+abs(fvaRange(:,2))==0);
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
    modelTemp.lb                   = fvaRange(:,1);
    modelTemp.lb(zeroRxns)         = [];
    modelTemp.ub                   = fvaRange(:,2);
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
model.prepTime    = elapsedTime;
model.rxnOrder = rxnOrder;
[~, model.rxnOrderRev] = sort(rxnOrder);
disp(['Done. Reactions removed ',num2str(length(zeroRxns))]);

%^ result structure
resInfo = struct();
resInfo.rxns = rxns;  %^ reactions tested
resInfo.rxnInLoops = rxnInLoops;  %^ reactions found to be in loops
resInfo.formulation = 0;  %^ formulation used (0 for Fast-SNP)
resInfo.nsTime = tSNP;  %^ Fast-SNP wall time
resInfo.nsCPU = cpuSNP;  %^ Fast-SNP CPU time
resInfo.loopPreprocessTime = 0;  %^ no loop preprocessing for Fast-SNP
resInfo.loopPreprocessCPU = 0;
resInfo.rxnLinkTime = 0;
resInfo.rxnLinkCPU = 0;
resInfo.lpTime = lpTime;  %^ wall time for solving LP problems (all zeros for Fast-SNP)
resInfo.milpTime = milpTime;  %^ wall time for solving MILP problems
resInfo.lpCPU = lpCPU;  %^ CPU time for solving LP problems (all zeros for Fast-SNP)
resInfo.milpCPU = milpCPU;  %^ CPU time for solving MILP problems
resInfo.nBinVar = nBinVar;  %^ number of binary variables in each solve
resInfo.mipStart = mipStart;  %^ whether MipStart is used in each solve
resInfo.optiSol = optiSol;  %^ whether the solution is optimal in each solve
resInfo.presolveOff = presolveOff;  %^ whether presolve is turned off in each solve
