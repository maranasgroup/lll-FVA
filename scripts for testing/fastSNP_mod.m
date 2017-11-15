function Nsnp = fastSNP_mod(model,solverName,params)
% Sparsify Nint solving a Null-space Basis Problem
% Inputs:     model structure, solverName (string, optional), params (solver parameters, optional)
% Outputs:    Nsnp (sparse null space matrix for internal loops)
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2016 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Structure parsing and matrix multiplication have been optimized by Joshua Chan
% (2017) for fair comparison with the method proposed along with localized
% loopless constraints
% (commented starting with '%^')
%%%%% Siu Hung Joshua Chan 2017

% Define solver parameters
if nargin<1
    disp('Not enough input parameters supported');return;
elseif nargin<2
    if exist('gurobi','file')  % Check available solvers in this order
        solverName = 'gurobi';
    elseif exist('cplexlp','file')
        solverName = 'cplex';
    else
        addpath('glpkmex');
        solverName = 'glpk';
    end
end
if nargin<3 && strcmp(solverName,'gurobi')
    params.outputflag     = 0;
    params.OptimalityTol  = 1e-9;
    params.FeasibilityTol = 1e-9;
elseif nargin<3 && strcmp(solverName,'cplex')
    options.Display       = 'off';
    options.TolFun        = 1e-9;
    options.TolRLPFun     = 1e-9;
elseif nargin<3 && strcmp(solverName,'glpk')
    params.tolobj         = 1e-9;
    params.tolbnd         = 1e-8;
    params.toldj          = 1e-8;
end

% Parse model and remove exchange rxns
model                     = parseInternalRxns_mod(model);
model.S(:,model.exchange) = [];
model.lb                  = -1e3*(model.lb(model.internal)<0);
model.ub                  = 1e3*(model.ub(model.internal)>0);

% Initialization of SNP parameters
% nullity = size(null(full(model.S)),2);  %^ This can be very time consuming
P_N     = 0;
Nsnp    = [];
fluxTol = 1e-6;
epsilon = 1e-3;

% Main loop
t0 = cputime;
disp('Performing Sparsest Null-space Persuit...')
weights = rand(1,size(model.S,2));   % Use a random nonnegative weight. Another suitable option is, weights = randn(1,size(model.S,2));

%^ create problem in advance and change in each loop to save time
modelSNP              = nullSparseBasisStructure_GUROBI_mod(model,zeros(1, size(model.S, 2)));

for jx = 1:min(size(model.S)) %^ nullity  % Iterate until we reached the desired number of basis vectors
    Ntemp    = [];
    
    %^ the multiplication below is very time consuming and unnecessary
    %^     P_NT     = speye(size(model.S,2))-P_N'*P_N;  
    %^     xJ = weights*P_NT;
    xJ = weights - (weights * P_N') * P_N;
    
    % (1) Run first LP optimization
    if strcmp(solverName,'gurobi')
        %^ No need to regenerate the same structure in each loop.
        %^         modelSNP              = nullSparseBasisStructure_GUROBI(model,weights*P_NT);
        %^ Just changing the last row is sufficient
        modelSNP.A(end, 1:size(model.S, 2)) = xJ;
        modelSNP.rhs(end)     = epsilon;
        modelSNP.sense(end)   = '>';
        sol                   = gurobi(modelSNP,params);
    elseif strcmp(solverName,'cplex')
        modelSNP              = nullSparseBasisStructure_CPLEX(model);
        modelSNP.options      = options;
        modelSNP.Aineq(end,:) = [-weights*P_NT,zeros(1,length(model.lb))];
        modelSNP.bineq(end)   = -epsilon;
        [xopt,~,exitFlag]     = cplexlp(modelSNP);
    else
        modelSNP              = nullSparseBasisStructure_GLPK(model,weights*P_NT);
        modelSNP.b(end)       = epsilon;
        modelSNP.ctype(end)   = 'G';
        [xopt,~,exitFlag]     = glpk(modelSNP.c,modelSNP.A,modelSNP.b,modelSNP.lb,modelSNP.ub,modelSNP.ctype,modelSNP.vartype,modelSNP.sense,params);
    end
    
    % Save solution (if possible)
    if (strcmp(solverName,'gurobi')&&strcmp(sol.status,'OPTIMAL'))||(strcmp(solverName,'cplex')&&(exitFlag==1))||(strcmp(solverName,'glpk')&&(exitFlag==5))
        if strcmp(solverName,'gurobi')
            vopt = sol.x(model.internal);
        else
            vopt = xopt(model.internal);
        end
        vopt(abs(vopt)<fluxTol) = 0;         % Set to zero entries below tol
        min_val = min(abs(vopt(vopt~=0)));
        Ntemp   = [Ntemp,vopt/min_val];
    end
    
    % (2) Run second LP optimization
    if strcmp(solverName,'gurobi')
        %^         modelSNP              = nullSparseBasisStructure_GUROBI(model,weights*P_NT);
        modelSNP.A(end, 1:size(model.S, 2)) = xJ;
        modelSNP.rhs(end)     = -epsilon;
        modelSNP.sense(end)   = '<';
        sol                   = gurobi(modelSNP,params);
    elseif strcmp(solverName,'cplex')
        modelSNP              = nullSparseBasisStructure_CPLEX(model);
        modelSNP.options      = options;
        modelSNP.Aineq(end,:) = [weights*P_NT,zeros(1,length(model.lb))];
        modelSNP.bineq(end)   = -epsilon;
        [xopt,~,exitFlag]     = cplexlp(modelSNP);
    else
        modelSNP              = nullSparseBasisStructure_GLPK(model,-weights*P_NT);
        modelSNP.b(end)       = epsilon;
        modelSNP.ctype(end)   = 'L';
        [xopt,~,exitFlag]     = glpk(modelSNP.c,modelSNP.A,modelSNP.b,modelSNP.lb,modelSNP.ub,modelSNP.ctype,modelSNP.vartype,modelSNP.sense,params);
    end
    
    % Save solution (if possible)
    if (strcmp(solverName,'gurobi')&&strcmp(sol.status,'OPTIMAL'))||(strcmp(solverName,'cplex')&&(exitFlag==1))||(strcmp(solverName,'glpk')&&(exitFlag==5))
        if strcmp(solverName,'gurobi')
            vopt = sol.x(model.internal);
        else
            vopt = xopt(model.internal);
        end
        vopt(abs(vopt)<fluxTol) = 0;         % Set to zero entries below tol
        min_val = min(abs(vopt(vopt~=0)));
        Ntemp   = [Ntemp,vopt/min_val];
    end
    
    % Select best solution
    if isempty(Ntemp)
        disp(['Premature ending. Some basis are unfeasible. Found: ',...
            num2str(jx-1), ' loop laws.']);
            %num2str(jx-1),'/',num2str(nullity),' loop laws.']);
        break;
    end
    [~,ix] = max(sum(Ntemp==0));       % Keep sparsest vector
    Nsnp = [Nsnp,Ntemp(:,ix)];
    
    % Re-build projetion onto col(x)
    P_N = orth(Nsnp)';
end
disp(['Pre-processing time: ',num2str(cputime-t0), ' cputime (s)']);