function [ initConditionForXMatrix ] = initConditionForxMatrix( n )
    tmpMatrix = eye(n, n);
    initConditionForXMatrix = reshape(tmpMatrix, [n*n, 1]);
end

