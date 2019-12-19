function [Zr, t] = fomp22(Y, A, A_t, C, N, k)
%% 2D orthogonal matching pursuit
% https://iataem.nwsuaf.edu.cn/yfang/
% [1]FANG Y, WU J, HUANG B. 2D sparse signal recovery via 2D orthogonal matching pursuit[J]. Science China Information Sciences, 2012, 55(4): 889897.
% [2]FANG Y. Sparse Matrix Recovery from Random Samples via 2D Orthogonal Matching Pursuit[J]. CoRR, 2011, abs/1101.5755.
% Zr = omp2(Y, A, A_t, C, N, k)
% z: vectorized recovered spikes
% Y: (m x m) sample matrix
% A: (m x n) sampling matrix, A = Phi * Psi
% A_t: A'
% C: (n x n) matrix for correlations between columns of A, i.e. C(i,j) = A_t(i,:) * A(:,j)
% N: (n x n) matrix for atom norms, i.e. N(i,j) = (C(i,i)*C(j,j))^0.5
% k: sparsity level

n = size(A,2);
Z0 = A_t * Y * A;       % 
R = Y;                  % residue matrix
F = ones(n);            % flag matrix
Lambda = []; % coordinates of selected atom
g = [];      %    

for t=1:k
    % find the most significant atom and record its coordinates
    P = abs(A_t * R * A ./ N) .* F;  % projection
    [x,I] = max(P); [x,j] = max(x);
    i = I(j);
    Lambda(t,1:2) = [i,j];    % record the coordinates
    F(i,j) = 0;             % clear the flag of the (i,j)-th atom
    % construct H and g, and calculate the optimal vector z
    if t==1
        A_11 = C(Lambda(1,1), Lambda(1,1)) .* C(Lambda(1,2), Lambda(1,2));
        A_11_inv = 1/A_11;
    else
        A_12 = C(Lambda((1:(t-1)),1), Lambda(t,1)) .* C(Lambda((1:(t-1)),2), Lambda(t,2));
        A_22 = C(Lambda(t,1), Lambda(t,1)) .* C(Lambda(t,2), Lambda(t,2));
        temp = (A_11_inv*A_12);                     % O((t-1)^2)
        F_22 = A_22 - ((A_12')*temp);               % O(t-1)
        F_11_inv = A_11_inv + temp*(temp')/F_22;    % O((t-1)^2)
        A_12_inv = -(F_11_inv*A_12)/A_22;           % O((t-1)^2)
        A_11_inv = [F_11_inv,  A_12_inv;
            A_12_inv', 1/F_22];
    end
    g(t,1) = Z0(Lambda(t,1), Lambda(t,2));
    z = A_11_inv * g(1:t);
    % update residue
    R = Y;
    for i=1:t
        R =  R - z(i) * A(:, Lambda(i,1)) * A_t(Lambda(i,2), :);
    end
end

Zr = sparse(Lambda(1:t,1), Lambda(1:t,2), z, n, n);
fprintf('%d spikes recovered.\n', t);