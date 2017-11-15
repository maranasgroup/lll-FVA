function modelMILP = looplessStructureMILP_GUROBI_mod(model,Nint)
% Builds loopless structure problem (gurobi)
% Inputs:  model structure, Nint
% Outputs: loopless model structure
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2016 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Edited the code to use sparse structure for faster parsing
% (commented starting with '%^')
%%%%% Siu Hung Joshua Chan 2017

[l,n] = size(model.S);

% Reorganize the stoichiometric matrix as [�nternalRxns|exchangeRxns]
modelMILP.numRxns  = n;
modelMILP.internal = model.internal;
modelMILP.exchange = model.exchange;
modelMILP.S        = model.S;
modelMILP.rxns     = model.rxns;
modelMILP.rxnNames = model.rxnNames;
modelMILP.mets     = model.mets;
modelMILP.metNames = model.metNames;
modelMILP.rev      = (model.lb<0).*(model.ub>0);
if ~isfield(model, 'description')
    model.description = '';
end
modelMILP.description = model.description;

% Calculate null basis of Sint (if necessary)
modelMILP.Nint = Nint';
[p,m]          = size(modelMILP.Nint);

% reactions in loops
rxnInLoops = find(any(modelMILP.Nint, 1));
% true number of binary variables
m = numel(rxnInLoops);

% Determine constrained internal reactions
M = sparse(1:m, rxnInLoops, 1, m, n);

% Define binary constant, model objective and sense
K                    = 1e3;
modelMILP.obj        = zeros(n+2*m,1);
modelMILP.sense      = 'max';

vUB = model.ub(modelMILP.internal(rxnInLoops));
vLB = model.lb(modelMILP.internal(rxnInLoops));

% Constraints definition
%^ use sparse matrix
modelMILP.A = [model.S,       sparse(l, 2 * m);...                                                 % Sv               = 0
               sparse(p, n),  modelMILP.Nint(:, rxnInLoops), sparse(p, m);...                      %    N'g           = 0
               sparse(m, n),  speye(m),                      sparse(1:m, 1:m, K + 1, m, m);...     %      g + (1+K)z <= K
               sparse(m, n), -speye(m),                      sparse(1:m, 1:m, -(K + 1), m, m);...  %     -g - (1+K)z <= -1
               M,             sparse(m, m), -sparse(1:m, 1:m, vUB, m, m); ...                      %  v     -   ub z <= 0
              -M,             sparse(m, m), -sparse(1:m, 1:m, vLB, m, m)];                         % -v     -   lb z <= -lb 
	
% RHS model
modelMILP.rhs = [zeros(l+p,1); K * ones(m, 1); -ones(m,1); zeros(m,1); -vLB];

milpInfo.con.Ng = (l + 1):(l + p);
milpInfo.con.gU = (l + p + 1):(l + p + m);
milpInfo.con.gL = (l + p + m + 1):(l + p + m * 2);
milpInfo.con.vU = (l + p + m * 2 + 1):(l + p + m * 3);
milpInfo.con.vL = (l + p + m * 3 + 1):(l + p + m * 4);

milpInfo.var.g = (n + 1):(n + m);
milpInfo.var.z = (n + m + 1):(n + m * 2);

% Sign assignation
%^ faster construction
modelMILP.sense = char(['=' * ones(1, l + p), '<' * ones(1, 4 * m)]);

% Assignation of variable types
%^ faster construction
modelMILP.vtype = char(['C' * ones(1, n + m), 'B' * ones(1, m)]);

% Bounds assignation 
modelMILP.lb  = [model.lb(1:n); -K * ones(m, 1); zeros(m,1)];
modelMILP.ub  = [model.ub(1:n);  K * ones(m, 1); ones(m,1)];

milpInfo.K = K;
modelMILP.milpInfo = milpInfo;
modelMILP.milpInfo.rxnLoopIds = zeros(n, 1);
modelMILP.milpInfo.rxnLoopIds(rxnInLoops) = 1:m;
[modelMILP.milpInfo.rxnFwdLoopIds, modelMILP.milpInfo.rxnRevLoopIds] = deal(modelMILP.milpInfo.rxnLoopIds);
