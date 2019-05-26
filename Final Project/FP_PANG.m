clc;clear;
close all
load S19_DTSSfit
%% Section 1: Controllability and observability
size(dtss_fit)

% Controllability
Q = ctrb(dtss_fit); % Controllability matrix
controllability = rank(Q)
% rank(Q) = N = 6 <--> Completely controllable

% Observability
R = obsv(dtss_fit);
observability = rank(R)
% rank(R) = N = 6 <--> Completely observable

%Check controllability for each single input
% for c = 1:2
%     sub_sys = dtss_fit(:,c);
%     ctrb_sub_sys(c) = rank(ctrb(sub_sys));
% end
% display(ctrb_sub_sys) %why the rank is 5
% 
% % Observability for each single output, and store the rank into an array.
% for c = 1:2
%     sub_sys = dtss_fit(c,:);
%     obsv_sub_sys(c) = rank(obsv(sub_sys));
% end
% display(obsv_sub_sys)
%% Section 2: LQR
q = 1;
% Q is N by N
Q = q * transpose(dtss_fit.C) * dtss_fit.C;

% R is M by M
R = 1E-1*diag([0.9 0.8]);
G = lqr(dtss_fit,Q,R);

Fs = 20;  % Sample rate (Hz)
Ts = 1/Fs; % Sample period (s)
T = 8; 
N = T / Ts; 
F = Fs / N; % frequency resolution
t = 0:Ts:T;

% for the first 1 sec
T1 = 1; % total period of data (s)
N1 = T1 / Ts; % total number of data

% Constant input excitation for the first 1s
u1 = 5;
u2 = -5;
x = []; y = []; u = [];
% u(:,1) = [u1; u2];
x(:,1) = zeros(6,1);
y(:,1) = zeros(1,2);

for k = 1:N1
    u(:,k) = [u1; u2];
    x(:,k + 1) = dtss_fit.A * x(:,k) + dtss_fit.B * u(:,k);
    y(:,k + 1) = dtss_fit.c * x(:,k) + dtss_fit.D * u(:,k);
end

% for the last 7 second
T2 = 8; % total period of data (s)
N2 = T2 / Ts; % total number of data

for k = (N1+1):N2
    u(:,k) = - G * x(:,k);
    x(:,k + 1) = dtss_fit.A * x(:,k) + dtss_fit.B * u(:,k);
    y(:,k + 1) = dtss_fit.c * x(:,k) + dtss_fit.D * u(:,k);
end
u(:,k + 1) = - G * x(:,k + 1);
%% Section 3
figure
subplot(2,1,1)
stairs(t,y')
line([0,t(end)],[0.5,0.5],'linestyle','--','color','m')
line([0,t(end)],[-0.5,-0.5],'linestyle','--','color','m')
grid on
ylabel('Output Signal (V)')
xlabel('Time (s)')
legend('y_1','y_2','0.5V','-0.5V')
subplot(2,1,2)
stairs(t,u')
grid on
ylabel('Input Signal (V)')
xlabel('Time (s)')
% ylim([-5.5 5.5])
legend('u_1','u_2')
set(gcf,'Color','white')
%% Section 4
% Process noise covariance matrix
alpha =5E-4; %process noise tuing gain 
QN = alpha*eye(2);

% Measurement noise covariance matrix
N = 72000; % Number of sample
u00= [zeros(1,N);zeros(1,N)];
yn = s19_plant(u00);
RN = cov(yn(1,:),yn(2,:));

% Cross covariance matrix between process noise and sensor noise
% Assume the sensor noise is uncorrelated with the process noise
NN = zeros(2,2);

[~,K] = kalman(dtss_fit,QN,RN,NN,'current'); 
disp('Kalman gain = ')
disp(K)
%% Section 5
tstart = 1.0;  tfinal = 8.0;  time = [tstart, tfinal]; 

% Output feedback control state space model
Aofc = [dtss_fit.A - K * dtss_fit.C];
Bofc = [dtss_fit.B - K * dtss_fit.D, K];
Cofc = -G;
dt_ofc = ss(Aofc,Bofc,Cofc,[]);

abs(eig(dt_ofc))

[y,u,xhat] = s19_plant(dt_ofc,time);
% y_hat = CofcT * xhat;
y_hat = dtss_fit.C * xhat + dtss_fit.D * u;

figure
subplot(3,1,1)
stairs(t,y(1,:))
hold on
stairs(t,y_hat(1,:))
stairs(t,abs(y_hat(1,:)-y(1,:)))
line([0,t(end)],[0.5,0.5],'linestyle','--','color','m')
line([0,t(end)],[-0.5,-0.5],'linestyle','--','color','m')
grid on
ylabel('Output Signal 2 (V)')
xlabel('Time (s)')
% legend('y_1','estimated y_1','0.5V','-0.5V')
legend({'$y_1$','$\hat{y}_1$','error'},'Interpreter','latex','FontSize',14,'Location','Eastoutside')
subplot(3,1,2)
stairs(t,y(2,:))
hold on
stairs(t,y_hat(2,:))
stairs(t,abs(y_hat(2,:)-y(2,:)))
line([0,t(end)],[0.5,0.5],'linestyle','--','color','m')
line([0,t(end)],[-0.5,-0.5],'linestyle','--','color','m')
grid on
ylabel('Output Signal 2 (V)')
xlabel('Time (s)')
legend({'$y_2$','$\hat{y}_2$','error'},'Interpreter','latex','FontSize',14,'Location','Eastoutside')
subplot(3,1,3)
% stairs(t,u')
% hold on
stairs(t,u')
ylabel('Input Signal (V)')
xlabel('Time (s)')
grid on
set(gcf,'Color','white')
legend({'$u_2$','$u_2$'},'Interpreter','latex','FontSize',14,'Location','Eastoutside')