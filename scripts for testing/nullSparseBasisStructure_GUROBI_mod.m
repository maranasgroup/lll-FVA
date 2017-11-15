function modelSNP = nullSparseBasisStructure_GUROBI_mod(modelSNP,x)
% Defines model structure for Fast-SNP
% Inputs:     model structure
%
% Outputs:    model structure for Fast-SNP (gurobi)
%%%%%%%%%%%%%%%%%%%%%% Pedro Saa UQ 2016 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Modified for faster parsing
% (commented starting with '%^')
%%%%% Siu Hung Joshua Chan 2017

% Parameter initialization
[m,n] = size(modelSNP.S);

% Definition of model structure
%^ use sparse matrix
modelSNP.A = [modelSNP.S, sparse(m,n);... 
    -speye(n), speye(n);...
    speye(n),speye(n);...
    x,sparse(1,n)];

% Redefine LP problem:
% (1) RHS
modelSNP.rhs = zeros(m+2*n+1,1);

% (2) Sense
%^ not using for loop
modelSNP.sense = char(['=' * ones(1, m), '>' * ones(1, 2*n + 1)]);
%^ modelSNP.sense = blanks(m+2*n+1);
%^ for ix = 1:m+2*n+1
%^     if ix<=m
%^         modelSNP.sense(ix) = '=';
%^     else
%^         modelSNP.sense(ix) = '>';
%^     end
%^ end

% (3) Bounds
modelSNP.lb = [modelSNP.lb;-1e3*ones(n,1)];
modelSNP.ub = [modelSNP.ub;1e3*ones(n,1)];

% (4) Model sense
modelSNP.modelsense = 'min';

% (5) Variables
modelSNP.vtype = 'C';

% (6) Obj fxn
modelSNP.obj = [zeros(n,1);ones(n,1)];