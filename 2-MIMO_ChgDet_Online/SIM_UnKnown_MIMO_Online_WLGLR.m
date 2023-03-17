clear;
clc;
tic
close all

%% Parameters
MC_Times = 1E7; % Monte Carlo simulation times

% 0- Scenario
Tx_PW_dBm = 23; % transmit power
AWGN_PW_dBmHz = -169; % AWGN power
BW_Hz = 1E7; % bandwidth
d_m = 100; % distance between UE and BS

% 1- enable signal
en_SIM = 1; % change or not change?  0: no change; 1: changed
en_EstMethod = 2; % estimation method  1: Shrinkage; 2:ML
KAPPA = 8; % condition number upper bound
xi = 40; % window size
xi_bar = 20;
theta = 120; % threshold

% 2- System Parameters
T = 8; % training length
Mt = 8; % number of TX antennas
Mr = 2; % number of RX antennas
N = 1E7; % frame length
SNR_dB = Function_Online.fun_generateSNR(Tx_PW_dBm,AWGN_PW_dBmHz,BW_Hz,d_m); % range of the transmit signal-to-noise ratio
CarrFreq_GHz = 80;

% 2-1 System Parameter(Trans unit of measurement)
SNR_1  = 10^(SNR_dB/10);
CarrFreq_Hz=CarrFreq_GHz * 1e9;
M = Mt*Mr;
WindowLength = xi - xi_bar;

% 3- Channel Parameters
PSI_degree = 30; % multipath parameters
PSI_pi = PSI_degree /360*2*pi;
dPer_timesLamda = 3;
delta_OMEGA_degree = 0.5; % describes how big the covariance matrix changes

% 4- Matrix Used in Channel Cov Generation
D_t = Function_Online.fun_generateDt(Mt,Mr);
D_r = Function_Online.fun_generateDr(Mr,Mt);

%% Main part
% Pilot
PILOT_MTX = Function_Online.fun_generatePilots(Mt,T);

% Noise PW
SIGMA_New = 1/SNR_1/T;

% generate channel
OMEGA_degree_1 = 0;
OMEGA_degree_2 = OMEGA_degree_1 + delta_OMEGA_degree;
[H_0,C_0] = Function_Online.fun_generateCovH(Mt,Mr,PSI_pi,CarrFreq_Hz,OMEGA_degree_1,dPer_timesLamda,D_t,D_r); % Channel covariance MTX in the previous frame
[H_1,C_1] = Function_Online.fun_generateCovH(Mt,Mr,PSI_pi,CarrFreq_Hz,OMEGA_degree_2,dPer_timesLamda,D_t,D_r); % Channel covariance MTX in the current frame (if changed)
C0_MTX = C_0 + SIGMA_New*eye(M);
C0_det = abs(det(C0_MTX));
C1_MTX = C_1 + SIGMA_New*eye(M);
C1_det = abs(det(C1_MTX));

% change point
if en_SIM == 0
    ChgPoint = inf;
else
    ChgPoint = xi;
end

% Initialization
Delay_set = ones(1,MC_Times)*-inf;
RunLen_set = ones(1,MC_Times)*-inf;
n = 1;
% Monte Carlo Simulation

parfor tT = 1:MC_Times
    en_change = 0;
    H_hat_set = zeros(M,xi);
    num = 1;
    DetChgFlag = 0;
    while num <= N
        Noise_MTX = Function_Online.fun_GenerateNoise(SNR_dB,Mr,T);
        
        % change time
        if num == ChgPoint
            en_change = 1;
        end
        
        % get channel
        if en_change == 0 % not changed
            H_0 = Function_Online.fun_GenerateH(C_0,Mt,Mr);
            Y = PILOT_MTX * H_0 + Noise_MTX;
        else                          % changed
            H_1 = Function_Online.fun_GenerateH(C_1,Mt,Mr);
            Y = PILOT_MTX * H_1 + Noise_MTX;
        end
        
        % parameter used in the estimation
        if en_SIM == 0
            SIGMA_B = Function_Online.fun_EstimateSIGMA_B(C_0);
        else
            SIGMA_B = Function_Online.fun_EstimateSIGMA_B(C_1);
        end
        
        % received signals
        h_hat = 1/T*PILOT_MTX'*Y;
        h_hat_vec = reshape(h_hat,[M,1]);
        if num < xi
            H_hat_set(:,num) = h_hat_vec;
        else
            H_hat_set(:,xi) = h_hat_vec;
        end
        
        % Whether have received enough data
        if num > xi
            w_num = 1;
            W_set = zeros(WindowLength,1);
            for w_Len = xi_bar:xi
                H_begin = xi - w_Len + 1;
                H_end = xi;
                H_cons = H_hat_set(:,H_begin:H_end);
                switch en_EstMethod
                    case 1
                        C_est = Function_Online.EstCovMTX_Shrinkage(H_cons,w_Len,T,SNR_dB,M); % Shrinkage-based estimation algorithm
                    case 2
                        C_est = Function_Online.EstCovMTX_ML(H_cons,w_Len,T,SNR_dB,SIGMA_B,KAPPA,M); % ML estimation algorithm
                end
                Cest_MTX = C_est + SIGMA_New*eye(M);
                Cest_det = abs(det(Cest_MTX));
                
                % LLR
                LLR0 = -w_Len*log(C0_det) - trace(C0_MTX\(H_cons*H_cons'));
                LLR1 = -w_Len*log(Cest_det) - trace(Cest_MTX\(H_cons*H_cons'));
                W = real(LLR1 - LLR0);
                
                % debug
                W_set(w_num) = W;
                w_num = w_num +1;
                
                if W > theta
                    DetChgFlag = 1;
                    break
                end
            end
            max_W = max(W_set);
            
            if DetChgFlag == 1
                break
            end
            
            % Update the H
            Mat = zeros(M,(xi-1));
            Mat = H_hat_set(:,2:xi);
            H_hat_set = zeros(M,xi);
            H_hat_set(:,1:(xi-1)) = Mat;
        end
        num = num + 1;
    end
    
    if en_SIM == 0 % false alarm
        RunLen_set(tT) = num - xi;
    else % detection delay
        Delay_set(tT) = num - ChgPoint;
    end
end

if en_SIM == 0 % false alarm
    ARL_set = mean(RunLen_set);
else % detection delay
    AverageDelay = mean(Delay_set);
end



toc
