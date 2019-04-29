function[bound, obasis, obfs, oval] = phaseTwo(A, b, c, sbasis, sbfs)
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
% sbasis a vector of size m of indices of column vectors for a feasible 
% basis for this problem from which to start the simplex method
% sbfs a vector of size n which is the basic feasible solution 
% corresponding to this basis
%
% Output:
% bound = 0 if the problem is unbounded (there is no optimal solution)
% bound = 1 if the problem is bounded (there is an optimal solution)
% obasis = a vector of size m of indices of column vectors which gives an 
% optimal feasible basis for the problem if the problem is bounded
% obfs = a vector of size n which is the optimal basic feasible solution 
% corresponding to this optimal basis if the problem is bounded
% oval = the objective value of this optimal basic feasible solution 
%(if the problem is bounded)

close all;
clc;

[m,n] = size(A);
obasis = zeros(m,1); %initializing return variables
obfs = zeros(n,1);
oval = 0;
nonbasic = createnonbasic(m,n,sbasis); %auxiliary vector to handle indices
 
AB = zeros(m,m);
CB = zeros(m,1);
AN = zeros(m,n-m);
CN = zeros (n-m,1);
r = zeros(n-m,1);
    

if m<n % Here, we make sure the matrix is not mxm to speed up process
    

    t = 0;
    while t == 0
        % t never changes from zero. the while breaks when 
                 % optimality or unboundedness are achieved
        

        for i = 1:m
            AB(:,i) = A(:,sbasis(i));
            CB(i) = c(sbasis(i));
        end

        for i = 1:n-m
            AN(:,i) = A(:,nonbasic(i));
            CN(i) = c(nonbasic(i));
        end
        
        Q = -AB\AN;
        p = AB\b;
        r = CN + (CB'*Q)'; 

        if sum(r(:)>0) == 0 %In this case, the optimum is achieved
            bound = 1;
            obasis = sbasis;
            for i = 1:m
                obfs(sbasis(i)) = p(i);
            end
            oval = c'*obfs;
            break; 
        end
        
        enter = 1; % We explore which variable enters the basis
        
        while   enter <= n-m && r(enter) <= 0
            enter = enter+1;
        end

        if sum(Q(:,enter)<0) == 0 % If this happens, the problem is unbounded
            bound = 0;
            oval = inf;
            obasis = sbasis;
            for i = 1:m
                obfs(sbasis(i)) = p(i);
            end
            break;
        end
        
        minimum = inf;
        leave = 1;
        for i = 1:m
            if Q(i,enter) < 0
                if -p(i)/Q(i,enter) < minimum
                    minimum =  -p(i)/Q(i,enter);
                    leave = i;
                end
            end
        end
        
        % Here we swap the values in the two sets. 
        sbasis(leave) = nonbasic(enter);
        % Now we reorder the entries
        sbasis = sort(sbasis);
        % Now we generate our new auxiliary vector
        nonbasic = createnonbasic(m,n,sbasis);
        
    end % we end while       
    
    
    
else % Here, if the matrix is mxm, we clear x and estimate objective value    
    obasis = (1:m)';
    x = A\b;
    bound = 1;
    oval = c'*x;  
    obfs = x;
    % Here, we don?t check for feasibility because phase I already did
    % that. The only basis possible here should be feasible allready.

end % we end if 


end % end function



function [nonbasic] = createnonbasic(m,n,sbasis)
nonbasic = zeros(n-m,1); %Vector of nonbasic indices
aux1 = 1; %Keeps index for basics while filling nonbasic
aux2 = 1; %Keeps index for nonbasics while filling nonbasic
%i is the current variable, it is either asigned to nonbasic or kept in  
%sbasis (do nothing)
for j=1:n
    if aux1 < m+1
        if sbasis(aux1) ~= j
            nonbasic(aux2) = j;
            aux2 = aux2+1;
        else
            aux1 = aux1+1;
        end
    else
        if aux2 < n-m+1
            nonbasic(aux2) = j;
            aux2 = aux2+1;
        end
    end
end
end
