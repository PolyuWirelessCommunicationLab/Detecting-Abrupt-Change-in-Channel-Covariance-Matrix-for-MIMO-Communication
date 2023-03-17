function [W_max, ChgPoint_est] = fun_Online_estW_UnknownC(en_EstMethod,H_hat_set, C0, C0_det,T,SNR_dB,KAPPA,SIGMA_B)
% H_hat_set : received hat_h;
% C0, C1: two covariance matries;
SNR_1  = 10^(SNR_dB/10);
SIGMA_New = 1/SNR_1/T;
[M,NUM] = size(H_hat_set);
W_set = zeros(1,NUM);

for n1 = 1:NUM
    H_ana = H_hat_set(:,n1:NUM);
    K = NUM-n1+1;
    
    % Est Covariance
    switch en_EstMethod
        case 1
            C_est = Function_Online.EstCovMTX_Shrinkage(H_ana,K,T,SNR_dB,M); % Shrinkage-based estimation algorithm
        case 2
            C_est = Function_Online.EstCovMTX_ML(H_ana,K,T,SNR_dB,SIGMA_B,KAPPA,M); % ML estimation algorithm
    end
    C1_MTX = C_est + SIGMA_New*eye(M);
    C1_det = abs(det(C1_MTX));
    
    LLR0 = -K*log(C0_det) - trace(C0\(H_ana*H_ana'));
    LLR1 = -K*log(C1_det) - trace(C1_MTX\(H_ana*H_ana'));
    W = (LLR1 - LLR0)/K;
    W_set(n1) = W;
end
[W_max, ChgPoint_est]  = max(W_set);
end