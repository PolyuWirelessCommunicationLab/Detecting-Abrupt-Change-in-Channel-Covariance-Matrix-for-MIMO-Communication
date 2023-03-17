clear;
clc;
tic
close all

%% Parameters
MC_Times = 1E6; % Monte Carlo simulation times

% 0- Scenario
Tx_PW_dBm = 23; % transmit power
AWGN_PW_dBmHz = -169; % AWGN power
BW_Hz = 1E7; % bandwidth
d_m = 20; % distance between UE and BS

% 1- enable signal
en_SIM = 1; % 0: no change; 1: changed
en_EstMethod = 1; % 1: Shrinkage; 2:ML
KAPPA = 60; % condition number upper bound
theta = -75;
LenSubInterval_1 = 10; % Length of sub-interval 1

% 2- System Parameters
T = 8; % training length
Mt = 8; % number of TX antennas
Mr = 64; % number of RX antennas
SNR_dB = Function_MIMO.fun_generateSNR(Tx_PW_dBm,AWGN_PW_dBmHz,BW_Hz,d_m);
CarrFreq_GHz = 80;

% 2-1 System Parameter(Trans unit of measurement)
SNR_1  = 10^(SNR_dB/10);
CarrFreq_Hz=CarrFreq_GHz * 1e9;
M = Mt*Mr;

% 3- Channel Parameters
PSI_degree = 50; % multipath parameters
PSI_pi = PSI_degree /360*2*pi;
dPer_timesLamda = 8;
delta_THETA_degree = 0.5; % describes how big the covariance matrix changes

% 4- Matrix Used in Channel Cov Generation
D_t = Function_MIMO.fun_generateDt(Mt,Mr);
D_r = Function_MIMO.fun_generateDr(Mr,Mt);

%% Main part
% Pilot
PILOT_MTX = Function_MIMO.fun_generatePilots(Mt,T);

% Noise PW
SIGMA_New = 1/SNR_1/T;

% Generate channel
load('C_1_Del0.mat') % load C_0
C_0 = C_1;
clear C_1
load(['C_1_Del',num2str(delta_THETA_degree),'.mat']) % load C_1
C0_MTX = C_0 + SIGMA_New*eye(M);
C0_det = abs(det(C0_MTX));

% parameter used in the estimation
SIGMA_B = 0.75;

% Monte Carlo Simulation
parfor tT = 1:MC_Times
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
        for det_l = 2 : LenSubInterval_1
            H_begin = LenSubInterval_1 - det_l + 1;
            H_end = LenSubInterval_1;
            K = det_l; % number of samples
            H_cons = H_hat_set(:, H_begin:H_end);
            
            % choose estimation method
            switch en_EstMethod
                case 1
                    C_est = Function_MIMO.EstCovMTX_Shrinkage(H_cons, K, T, SNR_dB, M); % Shrinkage-based estimation algorithm
                case 2
                    C_est = Function_MIMO.EstCovMTX_ML(H_cons,K,T,SNR_dB,SIGMA_B,KAPPA,M); % ML estimation algorithm
            end
            
            C1_est = C_est + SIGMA_New*eye(M);
            C1_det_est = abs(det(C1_est));
            
            LLR0 = -log(C0_det) - trace(C0_MTX\(H_cons*H_cons'));
            LLR1 = -log(C1_det_est) - trace(C1_est\(H_cons*H_cons'));
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
    RunLen_BigCirSet =  mean(RunLen_set);
else % detection delay
    Delay_BigCirSet = mean(Delay_set);
end

toc
