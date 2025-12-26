clc;
close all;
clear;
%%%%% Define functions for Claculating bandwidth & Envelope Detection%%%%% 
function[bandwidth] =calc_BW_baseband(M,f_vector) % 1% Criteria
M_magnitude = abs(M);
[maxi, ~] = max(M_magnitude);
threshold = 0.01 * maxi;
indices_above = find(M_magnitude >= threshold);
last_index = indices_above(end);
bandwidth = f_vector(last_index); 
end
function [BW] = calc_BW_passband(M, f_vector, fc)

M_mag = abs(M);

% Positive frequencies only
pos_idx = f_vector >= 0;
f_pos   = f_vector(pos_idx);
M_pos   = M_mag(pos_idx);

band_idx = (f_pos > 0) & (f_pos < fc);
f_band = f_pos(band_idx);
M_band = M_pos(band_idx);

% Peak inside sideband
peak_val = max(M_band);
threshold = 0.01 * peak_val;

idx = find(M_band >= threshold);

f_edge = f_band(idx(1));  
BW = fc - f_edge;
end

function [signaldemod] = env(signal, Vondiode, t, tau) % Envelope Detection
    N = length(signal);
    dt = t(2) - t(1); 
     v_cap = 0; % Initial voltage on the capacitor
       for i = 1:N
        if (signal(i) - Vondiode > v_cap)
            v_cap = signal(i) - Vondiode; 
        else
            v_cap = v_cap * exp(-dt/tau);
        end
         signaldemod(i) = v_cap;
      end
end
%-------------------%
%%%% Defining Parameters %%%%%%
fc = 50;
Ac = 1;
fm=18;  
ts=0.1/fc;  %Time step
T=10;      %Margin After Discontinuity at t=9;
fs=1/ts;     % Sampling frequency
N=ceil(T/ts); %Length of Vector
t=0:ts:(N-1)*ts;
df=fs/N;      %Frequency step
m=cos(2*pi*fm*t); %Message Signal
m(t>=9|t<=0)=0;
%--------------------% Plot Message in time
figure (1);
plot(t,m);
xlabel('Time (s)');
ylabel('m(t)');
title('Signal m(t) = cos(2*pi*18*t)');
axis tight;
grid on;
xlim([0 0.5])
%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Define Frequency vector
if(rem(N,2)==0)
f = - (0.5*fs) : df : (0.5*fs-df) ;
else
f = - (0.5*fs-0.5*df) : df : (0.5*fs-0.5*df) ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical fft
M_nu=fftshift(fft(m))*ts;

% Analytical expression
M_an = 4.5*sinc(9*(f-18)).*exp(-1i*pi*9*(f-18)) + 4.5*sinc(9*(f+18)).*exp(-1i*pi*9*(f+18));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot Message Spectrum
figure (2);
plot(f,abs(M_nu),"Color","b","LineWidth",1);
hold on 
plot(f,abs(M_an),"Color","r","LineWidth",1);
xlabel('Frequency (Hz)');
ylabel ('M(F)');
title(' Analytical and Numerical Fourier Transform');
axis tight;
grid on;
legend ("Numerical Fourier ", "Analytic FT")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%objective: calcukate bandwidth given a 1% criteria 
Band_Width=calc_BW_baseband(M_nu,f) 
%-------------------------------%
% Define Paramaters to Start Modulation 
c = Ac * cos(2*pi*fc*t); %carrier 
s =m.*c; % modulated signal 
S = fftshift(fft(s))*ts; % Spectrum of Modulated Signal

% Plot Modulated Signal Spectrum
figure(3)
plot(f,abs(S));
hold on ;
xlabel('Frequency (Hz)')
ylabel ('S(F)');
title('Modulated Signal ' )
grid on;
axis tight;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Using BPF : 1st method For SSB LSB
BPF = zeros(size(f));
BPF(f>(fc-30) & f<(fc))=1; %% +ve frequency
BPF(f<-(fc-30) & f>-(fc))=1;%% -ve frequency
S_LSB1 =  S.*BPF;
Band_Width_LSB1= calc_BW_passband(S_LSB1,f,fc)

%%%%%%%%%%%$$$$$$$$ 2nd method : hilbert filter 
m2=cos(2*pi*fm*t); % Message Signal
m2(t>=9|t<=0)=0;
mhilbert=imag(hilbert(m2)); 
c1=cos(2*pi*fc*t);
c2=sin(2*pi*fc*t);
ss =m2.*c1;
s_LSB2=ss+(mhilbert.*c2);
S_LSB2=fftshift(fft(s_LSB2))*ts;
Band_Width_LSB2= calc_BW_passband(S_LSB2,f,fc)

%----------------------------------%
%Plot Both LSB Methods
figure (4);
plot(f,abs(S_LSB1),'Color',"r");
hold on
plot(f,abs(S_LSB2),"Color","b" , "LineWidth",0.6);
xlabel('Frequency (Hz)');
ylabel ('|S(F)| LBS');
title('Modulated Message LSB 1st & 2nd Method');
axis tight;
grid on;
legend('1st Method: BPF', '2nd Method: Hilbert');

%----- Coherent Detector-----%%%%%%
s_LSB1 = ifft(fftshift(S_LSB1))/ts; % LSB in time domain
v = s_LSB1.*c;
V=fftshift(fft(v))*ts; % Demodulated Signal before LPF
H=abs(f)<fc; %LPF
low_pass = H.* V;
mo=real(ifft(fftshift(low_pass))/ts); % Demodulated Signal
figure(5)
plot(t,mo,"Color","b");
hold on
plot(t,m, "Color","r");
title('Signal & Demodulated SIgnal in Time Domain')
xlabel('Time');
ylabel('m(t)');
legend('Demodulated signal' , 'Orignal Message')
axis tight;
grid on; 
xlim([0 0.5])

%Amplitude Modulation Large Carrier
ka=0.4; % initial 
c =cos(2*pi*fc*t);
s_am =(1 + ka*m).*c; %Large Carrier Modulation                   
S_LC  =fftshift(fft(s_am))*ts;  % Spectrum
Tau=0.045;
g = env(s_am,0,t,Tau); %Demodulated signal using envelope detector
g_final= (g-mean(g))./ka;
%Large Carrier Spectrum
%%%% 6th Point (a)
figure(6)
plot(f, abs(S_LC), 'LineWidth', 0.7); 
hold on;
plot(f, abs(M_nu), 'LineWidth', 0.7);
grid on;
xlabel('frequency(hz)')
title('Large Carrier Modulation & Original Message : Spectrum')
legend('|S(f)|', '|M(f)|')
Band_Width_SLC= (calc_BW_passband(S_LC, f, fc)).*2

%%%% 6th Point (b)
figure(7)
plot(t, s_am, 'LineWidth', 0.7); 
hold on 
plot(t, g, 'LineWidth', 0.7);
xlabel('time(second)')
title('Envelope Detector Output & Original Message')
xlim([0 0.5])
legend('g(t)', '(1 + ka*m(t)*c)')
%%%% 6th Point (c)
figure(8)
plot(t, g_final, 'LineWidth', 0.7); 
hold on 
plot(t, m, 'LineWidth', 0.7);
grid on;
xlabel('time(second)')
title('Demodulated Signal & Orignal Message')
legend ("Demodulated Signal","Original Message")
xlim([0 0.5])

%%%% Figure 9 Spectrum 
G=fftshift(fft(g_final))*ts;
figure(9)
plot(f, abs(G), 'LineWidth', 0.7); 
hold on;
plot(f, abs(M_nu), 'LineWidth', 0.7);
grid on;
xlabel('Frequency (hz)')
title('Demodulated Signal & Original Message: Spectrum')
legend('G(f)', 'M(f)')

%%%%%%%%%%%5 6th Point (d)%%%%%%%%%%%%%%%%%%%%
ka=1.5; % overmodulated signal 
c =cos(2*pi*fc*t);
s =(1 + ka*m).*c;
s_am =(1 + ka*m).*c; %Large Carrier Modulation                   
S_LC  =fftshift(fft(s_am))*ts;  % Spectrum
Tau=0.045;
g = env(s_am,0,t,Tau);  %Envelope Output
g_final= (g-mean(g))./ka;   %Demodulated signal using envelope detector

%Large Carrier Spectrum
%6th Point (a) Revisited
figure(10)
plot(f, abs(S_LC), 'LineWidth', 0.7); 
hold on;
plot(f, abs(M_nu), 'LineWidth', 0.7);
grid on;
xlabel('frequency(hz)')
title('Large Carrier Modulation & Original Message : Spectrum')
legend('|S(f)|', '|M(f)|')
Band_Width_SLC2= (calc_BW_passband(S_LC, f, fc)).*2
%%%% 6th Point (b) revisited
figure(11)
plot(t, s_am, 'LineWidth', 0.7); 
hold on;
plot(t, g, 'LineWidth', 0.7);
grid on;
xlabel('time(second)')
title('Envelope Detector Output & Original Message')
xlim([0 0.5])
legend('g(t)', '(1 + ka*m(t)*c)')

