function x = getMipStart_gurobi(x0, N, conComp, zUB, options, K)
% presolve seems to be the problem
options.Presolve = 0;
if nargin < 3 || isempty(conComp)
    conComp = [ones(size(N, 1), 1); zeros(numel(x0) - size(N, 1), 1)];
end
if nargin < 4 || isempty(zUB)
    zUB = ones(size(N, 1), 1);
end
if nargin < 6
    K = 1000;
end
active = find(zUB);
rxnLoop = find(conComp > 0);
conCompActive = unique(conComp(rxnLoop(active)));
% faster than ismember this way
rxnActive = conComp(rxnLoop) == conCompActive(1);
for j = 2:numel(conCompActive)
    rxnActive = rxnActive | conComp == conCompActive(j);
end

n = numel(conComp);
if size(N, 1) < n
    N = [N; zeros(n - size(N, 1), size(N, 2))];
end

nLoop = numel(rxnLoop);
nAct = numel(active);
LP.A = [N(conComp > 0, :)',                      sparse(size(N, 2), nAct);  ...  % N'g = 0
        sparse(1:nAct, active, 1, nAct, nLoop),  (1 + K) * speye(nAct); ...  % g + (1 + K)z <= K
       -sparse(1:nAct, active, 1, nAct, nLoop), -(1 + K) * speye(nAct)];  % -g - (1 + K)z <= -1
LP.lb = [-1000 * rxnActive; zeros(nAct, 1)];
LP.ub = [1000 * rxnActive; ones(nAct, 1)];
LP.lb(x0(rxnLoop) < -1e-6 & zUB > 0) = 1;
LP.ub(x0(rxnLoop) > 1e-6 & zUB > 0) = -1;
LP.lb(nLoop + find(x0(rxnLoop(active)) > 1e-6)) = 1;
LP.ub(nLoop + find(x0(rxnLoop(active)) < -1e-6)) = 0;
LP.obj = zeros(nLoop + nAct, 1);
LP.rhs = [zeros(size(N, 2), 1); K * ones(nAct, 1); -ones(nAct, 1)];
LP.sense = char(['=' * ones(1, size(N, 2)), '<' * ones(1, nAct * 2)]);
LP.vtype = char(['C' * ones(1, nLoop), 'B' * ones(1, nAct)]);
sol = gurobi(LP, options);
z = zeros(nLoop, 1);
z(active) = sol.x((nLoop + 1) : end);
x = [x0; sol.x(1:nLoop); z];
end
