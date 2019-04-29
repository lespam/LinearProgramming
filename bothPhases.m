function[status, obasis, obfs, oval] = bothPhases(A, b, c)
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
% status = -1 if the feasible set is empty
% status = 0 if the feasible set is non-empty but the problem is unbounded 
%(there is no optimal solution)
% status = 1 if the problem is bounded (there is an optimal solution)
% obasis = a vector of size m of indices of an optimal feasible basis for 
% the problem if the feasible set is non-empty and the problem is bounded 
%(in terms of a set of indices of column vectors)
% obfs = a vector of size n which is the optimal basic feasible solution 
% corresponding to this optimal basis if the feasible set is non-empty and 
% the problem is bounded
% oval = the objective value of this optimal basic feasible solution (if 
% the feasible set is non-empty and the problem is bounded)

[nvac, basis, bfs] = phaseOne(A,b,c);
if nvac == 0
    status = -1;
    [m,n] = size(A);
    obasis = zeros(m,1);
    obfs = zeros(n,1);
    oval = 0;
else
    [status, obasis, obfs, oval] = phaseTwo(A, b, c, basis, bfs);
end

end