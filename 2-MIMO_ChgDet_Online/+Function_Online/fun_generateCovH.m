function [H,C] = fun_generateCovH(Mt,Mr,PSI_pi,CarrFreq,Theta_degree,dPer_timesLamda,Dt,Dr)
M = Mt*Mr;
%% generate C
lamda = 1/CarrFreq;
Theta_pi = Theta_degree/180*pi;
the_min = 0;
the_max = 2*pi;

dT = dPer_timesLamda*lamda*Dt;
dT_x = sin(Theta_pi)*dT;
dT_y = cos(Theta_pi)*dT;
dR = dPer_timesLamda*lamda*Dr;
dR_x = sin(Theta_pi)*dR;
dR_y = cos(Theta_pi)*dR;

fun = @(the) exp(-1i*2*pi/lamda*((1+PSI_pi^2/4*(-1+cos(2*the)))*dT_x + PSI_pi*sin(the)*dT_y + sin(the)*dR_x + cos(the)*dR_y));
C_1 = 1/2/pi*integral(fun,the_min,the_max,'ArrayValued',true);
C_2 = C_1';
C = (C_2+C_1)/2;

%% generate h
h_bo = sqrt(1/2)*(randn(M,1)+randn(M,1)*1j);
h = C^(1/2)*h_bo;
H = reshape(h,[Mt,Mr]);
