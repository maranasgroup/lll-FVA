function conComp = connectVars(N)
conComp = zeros(size(N, 1), 1);
conComp(~any(N, 2)) = -1;
nCon = 0;
vCur = false(size(N, 1), 1);
while ~all(conComp)
    vCur(:) = false;
    vCur(find(conComp == 0, 1)) = true;
    nCon = nCon + 1;
    nCur = 0;
    while nCur < sum(vCur)
        nCur = sum(vCur);
        vCur(any(N(:, any(N(vCur, :), 1)), 2)) = true;
    end
    conComp(vCur) = nCon;
end
