function s=SL0_2D(x, A, B, sigma_min )

    sigma_decrease_factor = 0.5;
    A_pinv = pinv(A);%О±Дж
    B_pinv = pinv(B);
    mu_0 = 2;
    L = 7;
s = A_pinv*x*(B_pinv)';
sigma = 2*max(abs(s));

% Main Loop
while sigma>sigma_min
    for i=1:L
        delta = OurDelta(s,sigma);
        s = s - mu_0*delta;
        s = s - A_pinv*(A*s*B'-x)*(B_pinv)';   % Projection
    end
    sigma = sigma * sigma_decrease_factor;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function delta=OurDelta(s,sigma)

delta = s.*exp(-s.^2/sigma.^2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SNR=estimate_SNR(estim_s,true_s)

err = true_s - estim_s;
SNR = 10*log10(sum(true_s.^2)/sum(err.^2));
end