function[nvac, basis, bfs] = phaseOne(A, b, c)
% @Leslie Pamela Brenes Valencia 143851
% @C?sar Becerra Campos 163256
%
% maximise c^T x
% subject to Ax = b, x >= 0, b >=0
% 
% Input:
% A mxn matrix with m <= n and rank of A is m
% b column vector with m rows
% c column vector with n rows
% 
% Output:
% nvac = 0 if the feasible set is empty
% nvac = 1 if the feasible set is non-empty
% basis = a vector of size m of indices of column vectors for a feasible 
% basis for the problem if the feasible set is non-empty
% bfs = a vector of size n of the basic feasible solution corresponding to 
% this basis (if the feasible set is non-empty)...

close all;
clc;

[m,n] = size(A);
basis = zeros(m,1);
bfs = zeros(n,1);


if m<n % Here, we make sure the matrix is not mxm to start finding a basis
    
    % We generate extended matrix and objective function
    B = zeros(m,m+n);
    C = zeros(m+n,1);
    B(1:m,1:n) = A;
    C = zeros(m+n,1); 
    C(n+1:n+m) = -1;
    
    %Here we shall use auxiliary variables (correction variables) that will
    %help us determine feasibility
    for i=1:m 
        if b(i) >= 0
            B(i,n+i) = 1;
        else
            B(i,n+i) = -1;
        end
    end
    
    sbasis = (n+1:n+m)';
    sbfs = zeros(m+n,1);
    sbfs(n+1:n+m) = abs(b);
    
    
    
    % We call phaseTwo to maximize the negative sum of the correction
    % variables, which maximizes when they are zero and they are all zero
    % if and only if the original problem is feasible
    
    [bound, obasis, obfs, oval] = phaseTwo(B,b,C,sbasis,sbfs);
    
    
    if oval < 0 %In this case the problem is not feasible
        nvac = 0;
    else % In this case, it is
        nvac = 1;
        basis = obasis;
        bfs = obfs(1:m);
    end
    
else % Here, since the matrix is mxm, we clear x and check if its feasible
    x = A\b;
    nvac = 1;
    for i = 1:n
        if x(i) < 0 % Here we check feasibility
            nvac = 0;
            basis (i) = i;
        end
    end
    if nvac == 1
        bfs = x;
        basis = [1:m]';
    end
end % we end if   

end % end function