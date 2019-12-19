function [Zr, t] = romp2d_v2_5(Y, A, A_t, C, N, k, ad)
%% 2D ROMP
% composed by Rinabell
% version 1.0 @18-04-28E
% version 2.0 @18-04-30
% version 2.5 @18-05-01

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
g = [];      
if nargin == 6
    ad = floor(n/10);
else if nargin < 6
        fprintf('input lack.\n');
    end
end

for t=1:k
    % find the most significant atom and record its coordinates
    P = abs(A_t * R * A ./ N) .* F;  % projection
    
    idx = sort(abs(P(:)),'descend');
    idx = idx(1:ad);
    for temp = 1:ad
        [i(temp),j(temp)] = find(P == idx(temp));
        F(i(temp),j(temp)) = 0;
        Lambda( (t-1) * ad + temp,1:2) = [i(temp),j(temp)];
    end
    
    if t==1
        A_11 = C(Lambda(1:(t*ad),1), Lambda(1:(t*ad),1)) .*...
            C(Lambda(1:(t*ad),2), Lambda(1:(t*ad),2));
        A_11_inv = pinv(A_11);
    else
        A_12 = C(Lambda((1: ((t-1)*ad) ), 1), Lambda( ((t-1)*ad+1):t*ad , 1)) .*...
            C(Lambda((1: ((t-1)*ad) ),2), Lambda( ((t-1)*ad + 1) :t*ad, 2));
        A_22 = C(Lambda( ((t-1)*ad+1):t*ad ,1), Lambda( ((t-1)*ad+1):t*ad ,1)) .*...
            C(Lambda( ((t-1)*ad+1):t*ad ,2), Lambda( ((t-1)*ad+1):t*ad ,2));
        temp = (A_11_inv*A_12);                     % O((t-1)^2)
        F_22 = A_22 - ((A_12')*temp);               % O(t-1)
        F_11_inv = A_11_inv + temp*pinv(F_22)*(temp');      %changed here    
        A_12_inv = -(F_11_inv*A_12)*pinv(A_22);     %changed here too
        A_11_inv = [F_11_inv,  A_12_inv; A_12_inv', pinv(F_22)];
    end
    for temp = ((t-1)*ad+1) :t*ad
        g(temp,1) = Z0( Lambda(temp,1),Lambda(temp,2) );
    end
    z = A_11_inv * g(1:t*ad);
    % update residue
    R = Y;
    for i=1:t*ad
        R =  R - z(i) * A(:, Lambda(i,1)) * A_t(Lambda(i,2), :);
    end
end
Zr = sparse(Lambda(1:t*ad,1), Lambda(1:t*ad,2), z, n, n);
fprintf('%d spikes recovered.\n', t);