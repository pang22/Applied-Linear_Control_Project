clc;clear;
close all
load S19_DTSSfit
%% LQR

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
%% Extra Credit

A = dtss_fit.A;
B = dtss_fit.B;
C = dtss_fit.C;
D = dtss_fit.D;
R = 1E-4*diag([0.199, 0.352]);

%initial states
x_hat_minus(:,1) = zeros(6,1);
x_hat(:,1) = zeros(6,1);
x_hat(:,2) = zeros(6,1);
% initial outputs
y(:,1) = zeros(2,1);

% initial error convariance
P_minus(:,:,1) = randn(6,6);

% u(:,1) = randn(2,1);
u(:,1) = zeros(2,1);

x(:,1) = zeros(6,1);
x(:,2) = zeros(6,1);
y = [];
for k = 2:N+1
    % nomial plant
    theta(:,k) = sqrt(R) * 0.1*randn(2,1);
    x(:,k) = A * x(:,k-1) + B * u(:,k-1);
    y(:,k) = C * x_hat(:,k-1) - D * u(:,k-1) + theta(:,k);
    
    % Prediction
    x_hat_minus(:,k) = A * x_hat_minus(:,k-1) + B * u(:,k-1);
    P_minus(:,:,k) = A * P_minus(:,:,k-1) * A'; % assume no processing noise
    
    % Correction
    K(:,:,k) = P_minus(:,:,k) * C' * inv(C * P_minus(:,:,k) * C' + R);
    x_hat = x_hat_minus + K(:,:,k) * (y(:,k) - C * x_hat_minus);
    P(:,:,k) = (eye(6) - K(:,:,k) * C) * P_minus(:,:,k);
    
    u(:,k) = - G * x_hat(:,k); % Dr. Southward told me that.
end

y_hat = C * x_hat + D * u;
%% Plot results

figure
set(gca,'fontsize',14) 
set(gcf,'Position',[0 0 800 600])
subplot(3,1,1)
stairs(t,y(1,:))
hold on
stairs(t,y_hat(1,:))
stairs(t,abs(y_hat(1,:)-y(1,:)))
grid on
ylabel('Output Signal 2 (V)')
xlabel('Time (s)')
legend({'$y_1$','$\hat{y}_1$','error'},'Interpreter','latex','FontSize',14,'Location','Eastoutside')
subplot(3,1,2)
stairs(t,y(2,:))
hold on
stairs(t,y_hat(2,:))
stairs(t,abs(y_hat(2,:)-y(2,:)))
grid on
ylabel('Output Signal 2 (V)')
xlabel('Time (s)')
legend({'$y_2$','$\hat{y}_2$','error'},'Interpreter','latex','FontSize',14,...
        'Location','Eastoutside')
subplot(3,1,3)
stairs(t,u')
ylabel('Input Signal (V)')
xlabel('Time (s)')
grid on
set(gcf,'Color','white')
legend({'$u_2$','$u_2$'},'Interpreter','latex','FontSize',14,...
       'Location','Eastoutside')

figure
set(gca,'fontsize',14) 
set(gcf,'Position',[0 0 800 600])
subplot(2,1,1)
hold on
grid on
for m = 1:2
    for p = 1:2
        plot(t,squeeze(K(p,m,:)))
    end
end
ylabel('Kalman Gain')
xlabel('Time (s)')
legend('Gain 1','Gain 2','Gain 3','Gain 4')
hold off

subplot(2,1,2)
hold on
grid on
for m = 1:6
    for p = 1:6
        plot(t,squeeze(P(p,m,:)))
    end
end
ylabel('P Matrix Element')
xlabel('Time (s)')
hold off