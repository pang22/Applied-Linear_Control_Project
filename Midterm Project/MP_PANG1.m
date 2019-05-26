%% Zhoubao Pang
%% Midterm Project
%% Applied Linear Control
%% Spring 2019
clc;clear;
close all
%% Section 1
Fs = 20; %Sample rate (Hz)
Fnyq = Fs/2;
Ts = 1/Fs;
N = 2^11; % Number of sample
F = Fs / N; % frequency resolution
T = N * Ts; %total period of data
t = Ts * [0:N-1];

% [b1,a1] = butter(2,[0.01 9]/(Fs/2)); %band pass filter
[b1,a1] = butter(2,[0.03 8]/(Fs/2)); %band pass filter

%linear chirp input
% u1 = 3*chirp(t,0.01,max(t),10,'linear');
% u2 = 3*chirp(t,0.01,max(t),10,'linear');

% filtered white noise input
u1 = filter(b1,a1,2.7*randn(1,N));
u2 = filter(b1,a1,2.7*randn(1,N));

u0 = zeros(1,N); %zero input

figure 
subplot(2,2,1)
nfft = 2^9;
wndo = nfft;
ovlp = nfft/2;
[Pxx,FR] = cpsd(u1,u1,wndo,ovlp,nfft,Fs);
semilogx(FR,10*log10(Pxx*Fs/2),'r')
xlim([0.01 Fs/2])
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
grid on
subplot(2,2,2)
plot(t,u1)
xlabel('Time (s)')
ylabel(' Input Voltage (v)')
grid on

subplot(2,2,3)
wndo = nfft;
ovlp = nfft/2;
[Pxx,FR] = cpsd(u2,u2,wndo,ovlp,nfft,Fs);
semilogx(FR,10*log10(Pxx*Fs/2),'r')
xlim([0.01 Fs/2])
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
grid on
subplot(2,2,4)
plot(t,u2)
xlabel('Time (s)')
ylabel(' Input Voltage (v)')
grid on
%% Section 2

u10= [u1;u0];
y1 = s19_plant(u10);

u01= [u0;u2];
y2 = s19_plant(u01);

% y1(1,:) = y10(1,:);
% y1(2,:) = y01(1,:);
% y2(1,:) = y10(2,:);
% y2(2,:) = y01(2,:);

u00= [u0;u0];
yn = s19_plant(u00);

SNR(1,1) = 10*log10((std(y1(1,:))^2 - std(yn(1,:))^2) / std(yn(1,:))^2); %dB
SNR(2,1) = 10*log10((std(y2(1,:))^2 - std(yn(1,:))^2) / std(yn(1,:))^2); %dB sensor 1 actuator 2
SNR(1,2) = 10*log10((std(y1(2,:))^2 - std(yn(2,:))^2) / std(yn(2,:))^2); %dB
SNR(2,2) = 10*log10((std(y2(2,:))^2 - std(yn(2,:))^2) / std(yn(2,:))^2); %dB
display(SNR)
%% Section 3
figure
[Pxxn1,FR] = pwelch(yn(1,:),wndo,ovlp,nfft,Fs);
[Pxxn2,~] = pwelch(yn(2,:),wndo,ovlp,nfft,Fs);

[Pxx11,~] = pwelch(y1(1,:),wndo,ovlp,nfft,Fs);
[Pxx12,~] = pwelch(y2(1,:),wndo,ovlp,nfft,Fs);
[Pxx21,~] = pwelch(y1(2,:),wndo,ovlp,nfft,Fs);
[Pxx22,~] = pwelch(y2(2,:),wndo,ovlp,nfft,Fs);

subplot(2,2,1)
semilogx(FR,10*log10(Pxx11*Fs/2))
hold on
semilogx(FR,10*log10(Pxx21*Fs/2))
semilogx(FR,10*log10(Pxxn1*Fs/2))
xlim([0.01 Fs/2])
xlabel('Frequency')
ylabel('Magnitude (dB)')
grid on
legend('Power11','Power21','n1','Location','southwest')
hold off

subplot(2,2,2)
semilogx(FR,10*log10(Pxx12*Fs/2))
hold on
semilogx(FR,10*log10(Pxx22*Fs/2))
semilogx(FR,10*log10(Pxxn2*Fs/2))
xlim([0.01 Fs/2])
xlabel('Frequency')
ylabel('Magnitude (dB)')
grid on
legend('Power12','Power22','n2','Location','southwest')
hold off

subplot(2,2,3)
plot(t,y1(1,:),t,y1(2,:),t,yn(1,:))
legend('y11','y21','n1','Location','southwest')
ylabel('Output Voltage (V)')
xlabel('Time (s)')
grid on

subplot(2,2,4)
plot(t,y2(1,:),t,y2(2,:),t,yn(2,:))
legend('y12','y22','n2','Location','southwest')
ylabel('Output Voltage (V)')
xlabel('Time (s)')
grid on
%% Section 4

% [Suu,~] = cpsd(u,u,wndo,ovlp,nfft,Fs);
%calculate Syy for weighting vector
[Syy(1,1,:),~] = cpsd(y1(1,:),y1(1,:),wndo,ovlp,nfft,Fs);
[Syy(1,2,:),~] = cpsd(y2(1,:),y2(1,:),wndo,ovlp,nfft,Fs);
[Syy(2,1,:),~] = cpsd(y1(2,:),y1(2,:),wndo,ovlp,nfft,Fs);
[Syy(2,2,:),~] = cpsd(y2(2,:),y2(2,:),wndo,ovlp,nfft,Fs);

% H1 = Suy./Suu;

[H1(1,1,:),fr] = tfestimate(u1,y1(1,:),wndo,ovlp,nfft,Fs);%first index is always sensor number
[H1(1,2,:),~] = tfestimate(u2,y2(1,:),wndo,ovlp,nfft,Fs);
[H1(2,1,:),~] = tfestimate(u1,y1(2,:),wndo,ovlp,nfft,Fs);
[H1(2,2,:),~] = tfestimate(u2,y2(2,:),wndo,ovlp,nfft,Fs);

gamma(1,1,:) = mscohere(u1,y1(1,:),wndo,ovlp,nfft,Fs);
gamma(1,2,:) = mscohere(u2,y2(1,:),wndo,ovlp,nfft,Fs);
gamma(2,1,:) = mscohere(u1,y1(2,:),wndo,ovlp,nfft,Fs);
gamma(2,2,:) = mscohere(u2,y2(2,:),wndo,ovlp,nfft,Fs);

% figure
% plot(FR,gamma1)
% hold on 
% plot(FR,gamma2)
% plot(FR,gamma3)
% plot(FR,gamma4)
% hold off
%% Section 5
% Nb = [3 4;...
%       2 4];
% Na = 4; %number of poles
% Nb = [4 5;...
%       5 5];
  
Nb = [4 5;...
      5 5];
Na = 5; %number of poles

% inte_s = tf([1], [1 0]); % CT integrator
% inte_z = c2d(inte_s,Ts,'Tustin'); %DT integrator
% b_inte = cell2mat(inte_z.Numerator);
% a_inte = cell2mat(inte_z.Denominator);
% Hinte = freqz(b_inte,a_inte,fr,Fs);
% Hinte = freqs([1],[1 0],2*pi*fr);
        
% weighting vector for invfreqz
wt = zeros(length(fr),1);
wt(2:50) = 1;
% wt(12:129) = 1;
%infreqz iteration
iter = 1E6;
        
for m = 1:2
    for p = 1:2
%         H1_wo_inte = squeeze(H1(p,m,:))./Hinte; % remove integerator
%         [num, den] = invfreqz(H1_wo_inte,(pi*fr)/(Fs/2),Nb(p,m),Na,[]);

%         wt = squeeze(Syy(p,m,:))/max(squeeze(Syy(p,m,:))); %use normalizedSyy as weighting vector
        
        [num, den] = invfreqz(squeeze(H1(p,m,:)),(pi*fr)/(Fs/2),Nb(p,m),Na,[],iter);

        %form cell array
        NUM{p,m,:} = num;
        DEN{p,m,:} = den;
    end
end

Htf = tf(NUM,DEN,Ts); 
% Htf = tf(NUM,DEN,Ts) * inte_z; %add integrator back
Hmr = minreal(Htf);
%% Section 6
% State-space realization
Hss = ss(Hmr); 

% Balanced realization
[Hbr, Hsv] = balreal(Hss); % [balance realization, hankel singular values]

% Number of states
NFIT = size(Hbr.A,1); % extract number of rows of A matrix
figure
bar([1:NFIT],Hsv)
grid on
set(gca,'yscale','log')
ylabel('Hankel Singular Value')
xlabel('State Index')

% State space model reduction
thresh = 0.1;  % based on visual inspection of the bar chart
elim = (Hsv < thresh);    % Boolean vector (true=eliminate) 
% elim = [zeros(5,1);ones(15,1)];
Hss_red = modred(Hbr, elim, 'Truncate'); 

max(abs(eig(Hss_red))) %test

%% Section 7
Poles = pole(Hss_red);
% Zeros = tzero(Hss_red);
[Wn, zeta] = damp(Hss_red);
Natural_Frequency_Hz = Wn / (2*pi);
T = table(Poles,Natural_Frequency_Hz,zeta);
disp(T)
% Section 8
figure
subplot(2,1,1)
plot(pole(Htf),'x','Color','r')
hold on
% plot(tzero(Htf),'o','Color','b')
hold off
legend('TF poles','AutoUpdate','off')
zgrid
axis equal
xlabel('RE (z)')
ylabel('IM (z)')

subplot(2,1,2)
plot(pole(Hss_red),'x','Color','r')
hold on
% plot(tzero(Hss_red),'o','Color','b')
hold off
legend('Reduced poles','AutoUpdate','off')
zgrid
axis equal
xlabel('RE (z)')
ylabel('IM (z)')
%% Section 9
FreqResp_tf = squeeze(freqresp(Htf,2*pi*fr)); 
Ph_tf = (180/pi)*angle(FreqResp_tf);
Mag_tf = 20*log10(abs(FreqResp_tf));

FreqResp_red = squeeze(freqresp(Hss_red,2*pi*fr)); 
Ph_red = (180/pi)*angle(FreqResp_red);
Mag_red = 20*log10(abs(FreqResp_red));
%% Section 10

% create 4-path plots
for m = 1:2
    for p = 1:2
        Name_str = sprintf('p=%d_m=%d',p,m);
        
        Ph_H1(p,m,:) = (180/pi)*angle(H1(p,m,:));
        Mag_H1(p,m,:) = 20*log10(abs(H1(p,m,:)));
        
        figure('Name',Name_str)
        subplot(3,1,1)
        semilogx(fr,squeeze(gamma(p,m,:)))
        grid on
        ylabel('Coherence')
        ylim([0 1])
        subplot(3,1,2)
        semilogx(fr,squeeze(Ph_H1(p,m,:)),fr,squeeze(Ph_tf(p,m,:)),fr,squeeze(Ph_red(p,m,:)),'-.');
        legend('H1','Htf','Hss red','Location','southwest')
        ylabel('Phase (deg)')
        ylim([-180 180])
        grid on
        subplot(3,1,3)
        semilogx(fr,squeeze(Mag_H1(p,m,:)),fr,squeeze(Mag_tf(p,m,:)),fr,squeeze(Mag_red(p,m,:)),'-.');
        ylabel('Mag (dB)')
        xlabel('Frequency (Hz)')
        legend('H1','Htf','Hss red','Location','southwest')
        grid on
    end
end
