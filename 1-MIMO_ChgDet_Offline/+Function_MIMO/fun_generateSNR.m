function SNR = fun_generateSNR(Tx_PW_dBm,AWGN_PW_dBmHz,BW_Hz,d_1)

d_km = d_1/1000;
L = -128.1 - 36.7*log10(d_km);
P_signal = Tx_PW_dBm + L;
P_AWGN = AWGN_PW_dBmHz + 10*log10(BW_Hz);
SNR = P_signal - P_AWGN;