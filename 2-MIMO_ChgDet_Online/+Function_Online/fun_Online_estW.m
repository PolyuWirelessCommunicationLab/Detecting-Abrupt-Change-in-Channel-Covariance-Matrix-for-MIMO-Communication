function [W_max, ChgPoint_est] = fun_Online_estW(H_hat_set, C0, C0_det, C1, C1_det)
% H_hat_set : received hat_h;
% C0, C1: two covariance matries;

[M,NUM] = size(H_hat_set);
Wsum_set = zeros(1,NUM);
W_sum = 0;
for n1 = NUM:-1:1
    H_ana = H_hat_set(:,n1);
    LLR0 = -log(C0_det) - trace(C0\(H_ana*H_ana'));
    LLR1 = -log(C1_det) - trace(C1\(H_ana*H_ana'));
    W = LLR1 - LLR0;
    W_sum = W_sum + W;
    Wsum_set(n1) = W_sum;
end
[W_max, ChgPoint_est]  = max(real(Wsum_set));
end