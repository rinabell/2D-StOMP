function [Zr, t] = stomp2d(Y, A, A_t, C, N, k)
%% 2D StOMP
% composed by Rinabell
% based on 2D-OMP by Fang Yong
% version 1.0 @18-05-03
% version 2.0E @18-05-06 with getting a cold
% version 3.0 @18-05-10
% version 4.0 @18-06-13
% version 4.5 @18-06-14
%% ini
% Zr = omp2(Y, A, A_t, C, N, k)
% z: vectorized recovered spikes
% Y: (m x m) sample matrix
% A: (m x n) sampling matrix, A = Phi * Psi
% A_t: A'
% C: (n x n) matrix for correlations between columns of A, i.e. C(i,j) = A_t(i,:) * A(:,j)
% N: (n x n) matrix for atom norms, i.e. N(i,j) = (C(i,i)*C(j,j))^0.5
% k: sparsity level (for m = 256, k =16 best)

[m,~] = size(A);
n = size(A,2);
Z0 = A_t * Y * A;       % 
R = Y;                  % residue matrix
F = ones(n);            % flag matrix
Lambda = []; % coordinates of selected atom
g = [];
if nargin <= 6
    k = floor(sqrt(m));
end
if nargin < 5
    fprintf('input lack.\n');
end
ts = 2;   % maybe try StOMP
list = zeros(k,1);

%% script
for t=1:k
    % find the most significant atom and record its coordinates
    P = abs(A_t * R * A ./ N) .* F;  % projection
    
    sigma = ts*norm(R)/sqrt(m);
    idx = sort(abs(P(:)),'descend');
    idxnum = max(find(idx>sigma));
    
    if t > 1
        if list(t-1) == 0 || (isempty(idxnum))
            break;
        end
    end
    
    if t==1
        
        list(t) = idxnum;
        for temp = 1:list(t)
            [i(temp),j(temp)] = find(P == idx(temp));
            F(i(temp),j(temp)) = 0;
            Lambda(temp,1:2) = [i(temp),j(temp)];
        end
        
        A_11 = C(Lambda(1:list(t),1), Lambda(1:list(t),1)) .*...
            C(Lambda(1:list(t),2), Lambda(1:list(t),2));
        A_11_inv = pinv(A_11);
        
    else
        
        list(t) = idxnum + list(t-1);
        for temp = 1:(list(t)-list(t-1))
            [i(temp),j(temp)] = find(P == idx(temp));
            F(i(temp),j(temp)) = 0;
            Lambda( list(t-1) + temp,1:2) = [i(temp),j(temp)];
        end
        
        A_12 = C(Lambda((1: list(t-1) ), 1), Lambda( (list(t-1)+1):list(t) , 1)) .*...
            C(Lambda((1: list(t-1)), 2), Lambda( (list(t-1) + 1) :list(t), 2));
        A_22 = C(Lambda((list(t-1)+1):list(t),1), Lambda((list(t-1)+1):list(t),1)) .*...
            C(Lambda((list(t-1)+1):list(t),2), Lambda((list(t-1)+1):list(t),2));
        temp = (A_11_inv*A_12);                     % O((t-1)^2)
        F_22 = A_22 - ((A_12')*temp);               % O(t-1)
        F_11_inv = A_11_inv + temp*pinv(F_22)*(temp');    
        A_12_inv = -(F_11_inv*A_12)*pinv(A_22);
        A_11_inv = [F_11_inv,  A_12_inv; A_12_inv', pinv(F_22)];
    end
    
    if t == 1
        for temp = 1:list(t)
            g(temp,1) = Z0( Lambda(temp,1),Lambda(temp,2) );
        end
    else
        for temp = (list(t-1)+1) : list(t)
            g(temp,1) = Z0( Lambda(temp,1),Lambda(temp,2) );
        end
    end
    z = A_11_inv * g(1:list(t));
    % update residue
    R = Y;
    for i=1:list(t)
        R =  R - z(i) * A(:, Lambda(i,1)) * A_t(Lambda(i,2), :);
    end
end

if t == 1
    Zr = sparse(Lambda(1:list(t-1),1), Lambda(1:list(t-1),2), z, n, n);
else
    Zr = sparse(Lambda(:,1), Lambda(:,2), z, n, n);
end
fprintf('%d spikes recovered.\n', t);