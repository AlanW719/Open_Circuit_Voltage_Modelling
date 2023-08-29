clear all; clc; close all;

data = xlsread('data.xls'); 
T = data(:,1); % Time in hours
I = data(:,2); % Currne in Amperes
V = data(:,3); % Voltage (V)
SOC = data(:,4); % SOC in [0 1]
N=length(SOC);

% SOC scaling
E=.175; % selected from lecture notes.
zs=SOC*(1-2*E)+E;
%%
%***********Question 2*************

%********Linear model*********
Po_L=[ones(N,1) SOC]; 
P_L = [Po_L I];
k_ls_L = (P_L'*P_L)\(P_L'*V); % LS estimate
R0est_L = k_ls_L(end); % Estimated ECM parameter R0
k_L = k_ls_L(1:2);
OCV_L = Po_L*k_L; % OCV of Linear model
V_L = P_L*k_ls_L; % Modeling voltage
e_L=V-V_L; % Voltage modeling error
% Output results
disp('OCV-SOC parameters of Linear Model:');
fprintf('a0 = %5.6f  a1 = %5.6f\n',k_L(1),k_L(2));
%%
% *********Polynomial model**********
Po_P(:,1)=ones(N,1);
Po_P(:,2)=zs;
iteration=0;
n_max=9; % Initial guess of the polynomial's highest order n 
m_max=8;  % Initial guess of the polynomial's highest order m 
index_n_m=zeros((n_max-1)*m_max,2);
R_sq_polynomial=zeros((n_max-1)*m_max,1);
St = sum((V-mean(V)).^2); % totoal sum of squares around the mean

% Test different order n and m in each iteration
for n=2:n_max
    for m=1:m_max
        iteration=iteration+1;
        Po_P(:,n+1)=zs.^n;
        Po_P(:,n+1+m)=zs.^(-m);
        PoP=Po_P(:,1:n+1+m); % remove extra column
        P_P=[PoP I];
        k_ls_P=pinv(P_P'*P_P)*P_P'*V; % LS estimate
        k_P=k_ls_P(1:length(k_ls_P)-1);
        OCV_P=PoP*k_P; % OCV of Polynomial model
        V_P = P_P*k_ls_P; % True voltage across battery terminals    
        Sr_polynomial = sum((V-V_P).^2);
        index_n_m(iteration,:)=[n m]; % index of the combination of order n and m in each iteration
        R_sq_polynomial(iteration) = (St-Sr_polynomial)/St;
    end
end

% Find optimum n and m of polynomial order
opt_idx = find(R_sq_polynomial==max(R_sq_polynomial)); % find the iteration value when the lowest RMSE occured
opt_n = index_n_m(opt_idx,1); % find the optimum order n
opt_m = index_n_m(opt_idx,2); % find the optimum order m

PoP_opt = PoP(:,[1:opt_n+1, n_max+2:n_max+1+opt_m]); % select the column relevant to the polynomial highest order opt_n and opt_m 
P_P_opt = [PoP_opt I];
k_ls_P_opt=pinv(P_P_opt'*P_P_opt)*(P_P_opt'*V); % LS estimate of the optimum polynomial
R0est_P = k_ls_P_opt(end); % Estimated ECM parameter R0
k_P_opt=k_ls_P_opt(1:length(k_ls_P_opt)-1);
OCV_P_opt=PoP_opt*k_P_opt; % OCV of optimum Polynomial model 
V_P_opt = P_P_opt*k_ls_P_opt; % True voltage across battery terminals
e_P=V-V_P_opt; % Modeling error (optimal order n and m have chosen for ploynomial model)
% Output results
fprintf('\n');
disp('OCV-SOC parameters of Polynomial Model:');
fprintf('optimal degree n = %d  optimal order m = %d\n',opt_n,opt_m);
for i=0:opt_n+opt_m
    fprintf('p%d = %5.6f\n',i,k_P_opt(i+1));
end
% Plot R square vs. iteration times
figure; box on;
x_i=1:1:iteration;
plot(x_i, R_sq_polynomial, 'LineWidth',1);
xlim([0 iteration]);
xlabel('Iteration'),ylabel('R^2');
title('R^2 of polynomial model in each iteration');
grid on;
%%
% *********Combined model**********
Po_C=[ones(N,1) 1./zs zs log(zs) log(1-zs)];
P_C=[Po_C I];
k_ls_C=(P_C'*P_C)\(P_C'*V); % LS estimate
R0est_C = k_ls_C(6); % Estimated ECM parameter R0
k_C=k_ls_C(1:5); % Estimated OCV parameters
OCV_C=k_C(1)*ones(length(zs),1) + k_C(2)*(1./zs) + ...
    k_C(3)*zs + k_C(4)*log(zs) + k_C(5)*log(1-zs); % OCV of Combined model
V_C = P_C*k_ls_C; % True voltage across battery terminals
e_C=V-V_C; % voltage error
% Output results
fprintf('\n');
disp('OCV-SOC parameters of Combined Model:');
for i=1:5
    fprintf('k%d = %5.6f\n',i-1,k_C(i));
end
%%
% **********Combined+3 model***********
Po_3=[ones(N,1) 1./zs 1./zs.^2 1./zs.^3 1./zs.^4 zs log(zs) log(1-zs)];
P_3=[Po_3 I];
k_ls_3=(P_3'*P_3)\(P_3'*V);
R0est_3 = k_ls_3(end); % Estimated ECM parameter R0
k_3=k_ls_3(1:8);
OCV_3=Po_3*k_3; % OCV of Combined+3 model
V_3 = P_3*k_ls_3; % True voltage across battery terminals
e_3=V-V_3; % voltage error
% Output results
fprintf('\n');
disp('OCV-SOC parameters of Combined+3 Model:');
for i=1:8
    fprintf('k%d = %5.6f\n',i-1,k_3(i));
end
%%
% *********Plot OCV vs. SOC*********
figure
plot(SOC,OCV_L, 'LineWidth', 1);
hold on;
plot(SOC,OCV_P_opt, 'LineWidth',1);
plot(SOC, OCV_C, 'LineWidth',1);
plot(SOC, OCV_3, 'LineWidth',1);
legend({'Linear model','Polynomial model', 'Combined model', 'Combined+3 model'}, 'Location', 'best');
xlim([0 1]);ylim([3.2 4.2]);xlabel('SOC'),ylabel('OCV (V)');
grid on;
%%
%***********Plot Modeling Error vs. SOC**************
idx_d = find(I<0); % discharging index
idx_c = find(I>0); % charging index
SOC_d = SOC(idx_d);
SOC_c = SOC(idx_c);
LinearModError_d = e_L(idx_d); % modeling error when discharging (Linear)
LinearModError_c = e_L(idx_c); % modeling error wehn charging (Linear)
PolynomialModError_d = e_P(idx_d); % modeling error when discharging (Polynomial)
PolynomialModError_c = e_P(idx_c); % modeling error wehn charging (Polynomial)
CombinedModError_d = e_C(idx_d); % modeling error when discharging (Combined)
CombinedModError_c = e_C(idx_c); % modeling error wehn charging (Combined)
Combined_3_ModError_d = e_3(idx_d); % modeling error when discharging (Combined+3)
Combiend_3_ModError_c = e_3(idx_c); % modeling error wehn charging (Combined+3)
%plot modeling error vs. SOC
figure;
subplot(221)
hold on; box on;
plot(SOC_d, LinearModError_d, 'LineWidth', 1);
plot(SOC_c, LinearModError_c, 'LineWidth', 1);
legend({'Discharge','Charge'}, 'Location', 'best');
xlim([0 1]);ylim([-0.35 0.1]);xlabel('SOC'),ylabel('Modeling Error (V)');
title('Linear Model');
grid on;

subplot(222)
hold on; box on;
plot(SOC_d, PolynomialModError_d, 'LineWidth', 1);
plot(SOC_c, PolynomialModError_c, 'LineWidth', 1);
legend({'Discharge','Charge'}, 'Location', 'best');
xlim([0 1]);ylim([-0.35 0.1]);xlabel('SOC'),ylabel('Modeling Error (V)');
title('Polynomial Model');
grid on;

subplot(223)
hold on; box on;
plot(SOC_d, CombinedModError_d, 'LineWidth', 1);
plot(SOC_c, CombinedModError_c, 'LineWidth', 1);
legend({'Discharge','Charge'}, 'Location', 'best');
xlim([0 1]);ylim([-0.35 0.1]);xlabel('SOC'),ylabel('Modeling Error (V)');
title('Combined Model');
grid on;

subplot(224)
hold on; box on;
plot(SOC_d, Combined_3_ModError_d, 'LineWidth', 1);
plot(SOC_c, Combiend_3_ModError_c, 'LineWidth', 1);
legend({'Discharge','Charge'}, 'Location', 'best');
xlim([0 1]);ylim([-0.35 0.1]);xlabel('SOC'),ylabel('Modeling Error (V)');
title('Combined+3 Model');
grid on;

% scaling (Referenced from lecuture)
s= 0:0.01:1;
Linear_error_c = zeros(1,length(s));
Linear_error_d = zeros(1,length(s));
Polynomial_error_c = zeros(1,length(s));
Polynomial_error_d = zeros(1,length(s));
Combined_error_c = zeros(1,length(s));
Combined_error_d = zeros(1,length(s));
Combined3_error_c = zeros(1,length(s));
Combined3_error_d = zeros(1,length(s));
for i = 2:length(s)
    id = find(SOC_d<s(i) & SOC_d>s(i-1)); % scaling index for discharging period
    ic = find(SOC_c<s(i) & SOC_c>s(i-1)); % scaling index for charging period
    Linear_error_d(i) = mean(LinearModError_d(id));
    Linear_error_c(i) = mean(LinearModError_c(ic));
    Polynomial_error_d(i) = mean(PolynomialModError_d(id));
    Polynomial_error_c(i) = mean(PolynomialModError_c(ic));
    Combined_error_d(i) = mean(CombinedModError_d(id));
    Combined_error_c(i) = mean(CombinedModError_c(ic));
    Combined3_error_d(i) = mean(Combined_3_ModError_d(id));
    Combined3_error_c(i) = mean(Combiend_3_ModError_c(ic));
end
Linear_error_both = (Linear_error_d.^2 + Linear_error_c.^2)/2;
Polynomial_error_both = (Polynomial_error_d.^2 + Polynomial_error_c.^2)/2;
Combined_error_both = (Combined_error_d.^2 + Combined_error_c.^2)/2;
Combined3_error_both = (Combined3_error_d.^2 + Combined3_error_c.^2)/2;

figure;
hold on; box on;
plot(s, Linear_error_both, 'LineWidth', 1);
plot(s, Polynomial_error_both, 'LineWidth', 1);
plot(s, Combined_error_both, 'LineWidth', 1);
plot(s, Combined3_error_both, 'LineWidth', 1);
legend({'Linear Model','Polynomial Model', 'Combined Model', 'Combined+3 Model'}, 'Location', 'best');
xlim([0 1]);ylim([0 0.065]);xlabel('SOC'),ylabel('Mean Squared Error (MSE)');
grid on;
%%
%***********Question 3*************
Sr_L = sum(e_L.^2); % Sum of fitting error squares (Sr) of Linear model
Sr_P = sum(e_P.^2); % Sum of fitting error squares (Sr) of Polynomial model
Sr_C = sum(e_C.^2); % Sum of fitting error squares (Sr) of Combined model
Sr_3 = sum(e_3.^2); % Sum of fitting error squares (Sr) of Combined+3 model

St = sum((V-mean(V)).^2); % totoal sum of squares around the mean
R_sq = (St-[Sr_L, Sr_P, Sr_C, Sr_3])/St; % colrelation coefficient of each of four models
figure
width = 0.4;
b = barh(R_sq, width, 'green');
xtips = b.YEndPoints+0.02;
ytips = b.XEndPoints;
labels1 = string(round(b.YData,4));
text(xtips,ytips,labels1,'VerticalAlignment','middle')
yticklabels({'Linear','Polynomial','Combined','Combined+3'});
xlabel('Correlation Coefficient R^2'); ylabel('Model'); title('Accuracy of each of four models');
xlim([0 1.16]);