clear all
close all
clc

%% Parameters
% 0- Scenario
Tx_PW_dBm = 23; % transmit power
AWGN_PW_dBm_per_Hz = -169; % AWGN power
BW_Hz = 1E7; % bandwidth
d_m = 100; % distance between UE and BS

% 2- System Parameters
T = 8; % training length
Mt = 8; % number of TX antennas
Mr = 2; % number of RX antennas
SNR_dB = Function_Online.fun_generateSNR(Tx_PW_dBm,AWGN_PW_dBm_per_Hz,BW_Hz,d_m); % range of the transmit signal-to-noise ratio
CarrFreq_GHz = 80; % carrier frequency

% 2-1-System Parameter (Trans unit of measurement)
SNR_1  = 10^(SNR_dB/10);
CarrFreq_Hz = CarrFreq_GHz * 1e9;
M = Mt*Mr;

% 3- Channel Covariance Parameters
PSI_degree = 30; % multipath parameters 
PSI_pi = PSI_degree /360*2*pi;
dPer_timesLamda = 3;
delta_OMEGA = 1; % describes how big the covariance matrix changes

% 4- Matrices Used in Channel Cov Generation
D_t = Function_Online.fun_generateDt(Mt,Mr);
D_r = Function_Online.fun_generateDr(Mr,Mt);


%% Main
% New Noise Power
SIGMA_New = 1/SNR_1/T;

% Generate Channel Covariance Matrix
OMEGA_degree_1 = 0;
OMEGA_degree_2 = OMEGA_degree_1 + delta_OMEGA;
[H_0,C_0] = Function_Online.fun_generateCovH(Mt,Mr,PSI_pi,CarrFreq_Hz,OMEGA_degree_1,dPer_timesLamda,D_t,D_r); % Channel covariance MTX in the previous frame
[H_1,C_1] = Function_Online.fun_generateCovH(Mt,Mr,PSI_pi,CarrFreq_Hz,OMEGA_degree_2,dPer_timesLamda,D_t,D_r); % Channel covariance MTX in the current frame (if changed)

C0_MTX = C_0 + SIGMA_New*eye(M);
C1_MTX = C_1 + SIGMA_New*eye(M);

% the analysis results of the Theorem 1
Phi_01 = real(-M - logdet(C1_MTX) +logdet(C0_MTX) + trace(C1_MTX/C0_MTX))


function y = logdet(A)
U = chol(A);
y = 2*sum(log(diag(U)));
end