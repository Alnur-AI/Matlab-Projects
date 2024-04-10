function [ dftx ] = findDftx( ftx )
%     dftx = zeros(0);
%     symArray = sym('x%d', [1 size(ftx, 1)]);
%     syms(symArray);
%     dftx = jacobian(ftx, symArray);
    dftx = zeros(0);
    n = size(ftx, 1);
    symArray = sym('x%d', [1 n]);
    syms(symArray);

    tmpArray = [];
    for i = 1:n
        tmpArray = [tmpArray, str2sym(ftx(i))];
    end
    dftx = jacobian(tmpArray, symArray);
end
