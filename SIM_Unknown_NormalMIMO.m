clear;
clc;
tic
close all

%% Parameters
MC_Times = 3; % Monte Carlo simulation times

% 0- Scenario
Tx_PW_dBm = 23;
AWGN_PW_dBmHz = -169;
BW_Hz = 1E7;
d_m = 20;

% 1- enable signal
en_SIM = 0; % 0: no change; 1: changed
en_EstimateMethod = 1; % 1: Shrinkage; 2:ML
KAPPA = 8; % condition number upper bound
theta = 10;
LenSubInterval_1 = 30; % Length of sub-interval 1

% 2- System Parameters
T = 8; % training length
Mt = 8; % number of TX antennas
Mr = 2; % number of RX antennas
SNR_dB = Function_MIMO.fun_generateSNR(Tx_PW_dBm,AWGN_PW_dBmHz,BW_Hz,d_m);
CarrFreq_GHz = 80;

% 3 System Parameter(Trans unit of measurement)
SNR_1  = 10^(SNR_dB/10);
CarrFreq_Hz=CarrFreq_GHz * 1e9;
M = Mt*Mr;

% 4- Channel Parameters
PSI_degree = 30; % multipath parameters —— Angle spread parameters (0.1,1,5) degree
PSI_pi = PSI_degree /360*2*pi;
dPer_timesLamda = 3;
delta_OMEGA_degree = 1; % How much channel changed

% 5- Matrix Used in Channel Cov Generation
D_t = Function_MIMO.fun_generateDt(Mt,Mr);
D_r = Function_MIMO.fun_generateDr(Mr,Mt);

% 6- Frame
Win_at_less = 15; % the smallest window length (to increase the estimation performance)


%% Main part
% Pilot
PILOT_MTX = Function_MIMO.fun_generatePilots(Mt,T);

% 3-2 Noise PW
SIGMA_New = 1/SNR_1/T;

% 3-3 Generate channel
OMEGA_degree_1 = 0;
OMEGA_degree_2 = OMEGA_degree_1 + delta_OMEGA_degree;
% Channel covariance MTX in the previous frame
[H_0,C_0] = Function_MIMO.fun_generateCovH(Mt,Mr,PSI_pi,CarrFreq_Hz,OMEGA_degree_1,dPer_timesLamda,D_t,D_r);
% Channel covariance MTX in the current frame (if changed)
[H_1,C_1] = Function_MIMO.fun_generateCovH(Mt,Mr,PSI_pi,CarrFreq_Hz,OMEGA_degree_2,dPer_timesLamda,D_t,D_r);
C0_MTX = C_0 + SIGMA_New*eye(M);
C0_det = abs(det(C0_MTX));

% parameter used in the estimation
if en_SIM == 0
    SIGMA_B = Function_MIMO.fun_EstimateSIGMA_B(C_0);
else
    SIGMA_B = Function_MIMO.fun_EstimateSIGMA_B(C_1);
end

% Initialization
Delay_set = zeros(1,MC_Times);
RunLen_set = zeros(1,MC_Times);

% Monte Carlo Simulation
for tT = 1:MC_Times
    % Generate the change point
    if en_SIM == 0
        ChgPoint = inf;
    else
        ChgPoint = randperm(LenSubInterval_1,1); % the change could happen at everywhere
    end
    
    % Generate the samples
    en_change = 0;
    H_hat_set = zeros(M,LenSubInterval_1);
    DetChgFlag = 0;
    for num =1 : LenSubInterval_1
        Noise_MTX = Function_MIMO.fun_GenerateNoise(SNR_dB,Mr,T);
        
        % change time
        if num == ChgPoint
            en_change = 1;
        end
        
        % get channel
        if en_change == 0
            % not changed
            H_0 = Function_MIMO.fun_GenerateH(C_0,Mt,Mr);
            Y = PILOT_MTX * H_0 + Noise_MTX;
        else
            % changed
            H_1 = Function_MIMO.fun_GenerateH(C_1,Mt,Mr);
            Y = PILOT_MTX * H_1 + Noise_MTX;
        end
        
        % received signals
        h_hat = 1/T*PILOT_MTX'*Y;
        h_hat_vec = reshape(h_hat,[M,1]);
        H_hat_set(:,num) = h_hat_vec;
        
    end
    
    % Detect Chg
    Num_OtherSubInt = 1;
    Flag_DetChg = 0;
    while Flag_DetChg == 0
        for det_l = (1+Win_at_less) : LenSubInterval_1
            H_begin = LenSubInterval_1 - det_l + 1;
            H_end = LenSubInterval_1;
            K = det_l; % number of samples
            H_cons = H_hat_set(:, H_begin:H_end);
            
            % choose estimation method
            switch en_EstimateMethod
                case 1
                    C_est = Function_MIMO.EstCovMTX_Shrinkage(H_cons, K, T, SNR_dB, M); % Shrinkage-based estimation algorithm
                case 2
                    C_est = Function_MIMO.EstCovMTX_ML(H_cons,K,T,SNR_dB,SIGMA_B,KAPPA,M); % ML estimation algorithm
            end
            
            C1_est = C_est + SIGMA_New*eye(M);
            C1_det_est = abs(det(C1_est));
            
            LLR0 = -det_l*log(C0_det) - trace(C0_MTX\(H_cons*H_cons'));
            LLR1 = -det_l*log(C1_det_est) - trace(C1_est\(H_cons*H_cons'));
            W = real(LLR1-LLR0);
            
            if W > theta
                Flag_DetChg = 1;
                break
            end
        end
        
        if Flag_DetChg == 0
            % generate the next covariance interval
            H_hat_set = zeros(M,LenSubInterval_1);
            for num =1 : LenSubInterval_1
                Noise_MTX = Function_MIMO.fun_GenerateNoise(SNR_dB,Mr,T);
                
                % get channel
                if en_change == 0
                    % not changed
                    H_0 = Function_MIMO.fun_GenerateH(C_0,Mt,Mr);
                    Y = PILOT_MTX * H_0 + Noise_MTX;
                else
                    % changed
                    H_1 = Function_MIMO.fun_GenerateH(C_1,Mt,Mr);
                    Y = PILOT_MTX * H_1 + Noise_MTX;
                end
                
                % received signals
                h_hat = 1/T*PILOT_MTX'*Y;
                h_hat_vec = reshape(h_hat,[M,1]);
                H_hat_set(:,num) = h_hat_vec;
                
            end
            Num_OtherSubInt = Num_OtherSubInt + 1;
        end
    end
    
    if en_SIM == 0 % false alarm
        RunLen_set(tT) = Num_OtherSubInt*LenSubInterval_1;
    else % detection delay
        Delay_set(tT) = Num_OtherSubInt*LenSubInterval_1 - ChgPoint;
    end
    
end

if en_SIM == 0 % false alarm
    ARL_mean =  mean(RunLen_set);
else % detection delay
    CADD_mean = mean(Delay_set);
end


toc
