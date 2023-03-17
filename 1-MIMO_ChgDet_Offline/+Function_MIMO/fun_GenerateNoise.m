function Noise_MTX = fun_GenerateNoise(SNR_dB,Mr,T)
SNR_1 = 10^(SNR_dB/10);
noise_PW = 1/SNR_1;
noise = sqrt(noise_PW/2)*(randn(Mr*T, 1)+1j*randn(Mr*T, 1));
Noise_MTX = reshape(noise,[T,Mr]);
end