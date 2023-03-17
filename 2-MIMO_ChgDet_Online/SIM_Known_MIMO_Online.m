
clear;
close all
clc;
tic

%% Parameters
MC_Times = 1E7; % Monte Carlo simulation times
theta = 1; % threshold

% 0- Scenario
Tx_PW_dBm = 23; % transmit power
AWGN_PW_dBmHz = -169; % AWGN power
BW_Hz = 1E7; % bandwidth
d_m = 100; % distance between UE and BS

% 1- enable signal
en_SIM = 0; % 0: no change; 1: changed

% 2- System Parameters
T = 8; % training length
Mt = 8; % number of TX antennas
Mr = 2; % number of RX antennas
N = inf; % frame length
SNR_dB = Function_Online.fun_generateSNR(Tx_PW_dBm,AWGN_PW_dBmHz,BW_Hz,d_m); % range of the transmit signal-to-noise ratio
CarrFreq_GHz = 80; % carrier frequency

% 2-1 System Parameter(Trans unit of measurement)
SNR_1  = 10^(SNR_dB/10);
CarrFreq_Hz=CarrFreq_GHz * 1e9;
M = Mt*Mr;

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
    ChgPoint = inf; % change not happens
else
    ChgPoint = 1; % change happens at the first coherence time interval
end

% Initialization
RunLen_set  = zeros(MC_Times,1);
Delay_set  = zeros(MC_Times,1);
parfor tT = 1:MC_Times
    en_change = 0; % to label whether the covariance matrix has changed in the considered interval (initialized at 0)
    DetChgFlag = 0;
    num = 1;
    W = 0;
    
    while num <= N
        % noise matrix
        Noise_MTX = Function_Online.fun_GenerateNoise(SNR_dB,Mr,T);
        
        % change time
        if num == ChgPoint
            en_change = 1; % if the interval is the ChgPoint, then the enable signal = 1
        end
        
        % get channel
        if en_change == 0 % not changed
            H_0 = Function_Online.fun_GenerateH(C_0,Mt,Mr);
            Y = PILOT_MTX * H_0 + Noise_MTX;
        else                     % changed
            H_1 = Function_Online.fun_GenerateH(C_1,Mt,Mr);
            Y = PILOT_MTX * H_1 + Noise_MTX;
        end
        % received signals
        h_hat = 1/T*PILOT_MTX'*Y;
        h_hat_vec = reshape(h_hat,[M,1]);
        
        LLR0 = -log(C0_det) - trace(C0_MTX\(h_hat_vec*h_hat_vec'));
        LLR1 = -log(C1_det) - trace(C1_MTX\(h_hat_vec*h_hat_vec'));
        W = real(LLR1-LLR0) + W; % sum of LLRs
        
        if W < 0
            W = 0;
        elseif W > theta
            break
        end
        num = num + 1;
        
    end
    
    if en_SIM == 0 % false alarm
        RunLen_set(tT) = num - 1;
    else               % detection delay
        Delay_set(tT) = num - 1;
    end
end

if en_SIM == 0 % false alarm
    ARL = mean(RunLen_set);
else                % detection delay
    CADD = mean(Delay_set);
end

toc
