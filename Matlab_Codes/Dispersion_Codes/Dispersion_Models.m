% eps(lambda) = 1 + sum (B_n * lambda^2)/(lambda^2 - C_n) %
% w_n = 1/sqrt(c_n) %
% gamma_n = 0 %
% sigma_n = B_n %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Sellmier model for BTO
% B = 4.187;
% C = 0.223;
% 
% w_0 = 1 ./ C; % Peak frequency
% gamma_n = 0;
% sigma_n = B;
% 
% View_Dispersion(w_0,gamma,sigma,1,0.4,0.7,500,"BTO")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Sellmier model for fused silica
% B = [0.6961663; 0.4079426; 0.8974794];
% C = [0.0684043; 0.1162414; 9.896161];
% 
% w_0 = 1 ./ C; % Peak frequency
% gamma_n = [0; 0; 0];
% sigma_n = B;
% eps_infin = 1;
% 
% View_Dispersion(w_0,gamma_n,sigma_n,eps_infin,0.21,6.7,1500,"SiO2")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Meep model for SiO2
w_0 = 1 ./ 0.103320160833333; % Peak frequency
gamma_n = 1 ./ 12.3984193000000;
sigma_n = 1.12;
eps_infin = 1;

View_Dispersion(w_0,gamma_n,sigma_n,eps_infin,0.25,1.77,800,"SiO2")
