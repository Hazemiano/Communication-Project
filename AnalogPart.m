clc;
%pkg load signal
close all;
clear;
function[bandwidth] =calc_BW(M,f_vector)
M_magnitude = abs(M);
[maxi, index] = max(M_magnitude);
threshold = 0.01 * maxi;
indices_above = find(M_magnitude >= threshold);
last_index = indices_above(end);
bandwidth = f_vector(last_index); 
end
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

figure (1);
plot(t,m);
xlabel('Time (s)');
ylabel('m(t)');
title('Signal m(t) = cos(2*pi*18*t)');
axis tight;
grid on;
%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%plot
figure (2);
plot(f,abs(M_nu),"Color","b","LineWidth",1);
xlabel('Frequency (Hz)');
ylabel ('M(F)');
title(' Analytical and Numerical Fourier Transform');
axis tight;
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%objective: calcukate bandwidth given a 1% criteria 
Band_Width=calc_BW(M_nu,f)
% The rest of this code should clculate the bandwidth 1% criteria 
 %graphically BW=21.6hz , hint: check zerocrossings function 
 c = Ac * cos(2*pi*fc*t); 
 s =m.*c;
S = fftshift(fft(s))*ts;
figure(4)
plot(f,abs(S));
hold on ;
xlabel('Frequency (Hz)')
ylabel ('S(F)');
title('Modulated Signal ' )
grid on;
axis tight;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Using BPF : 1st method 
BPF = zeros(size(f));
BPF(f>(fc-Band_Width) & f<(fc))=1; %% +ve frequency
BPF(f<-(fc-Band_Width) & f>-(fc))=1;%% -ve frequency
S_LSB1 =  S.*BPF;
Band_Width_LSB1= calc_BW(S_LSB1,f)
%%%%%%%%%%%$$$$$$$$ 2nd method : Ignore Lpf , use hilbert filter 
LPF = abs(f) <fc; %%%%LPF
S_LSB2 = LPF.* S;
Band_Width_LSB2= calc_BW(S_LSB2,f)

figure (5);
plot(f,abs(S_LSB1),'Color',"r");
hold on
plot(f,abs(S_LSB2),"Color","b" , "LineWidth",0.6);
xlabel('Frequency (Hz)');
ylabel ('|S(F)| LBS');
title('Modulated Message 1st & 2nd Method');
axis tight;
grid on;
legend('1st Method LSB', '2nd Method LSB');

s_LSB1 = ifft(fftshift(S_LSB1))/ts; % LSB in time domain
v = s_LSB1.*c;
V=fftshift(fft(v))*ts; % Demodulated Signal before LPF
H=abs(f)<fc; %LPF
low_pass = H.* V;
figure(6)
plot(f,abs(low_pass))
grid on
mo=real(ifft(fftshift(low_pass))/ts); % Demodulated Signal
figure(7)
plot(t,mo);
hold on
plot(t,m, "Color","r");
title('Signal & Demodulated SIgnal in Time Domain')
xlabel('Time');
ylabel('m(t)');
legend('Demodulated signal' , 'Orignal Message')
axis tight;
grid on;
%Amplitude Modulation Large Carrier
ka=0.4; % initial 
c =cos(2*pi*fc*t);
s =(1 + ka*m).*c;                
g  =abs(hilbert(s)); %envelope detector output (demodulated signal)        
S_LC  =fftshift(fft(s))*ts;   
am =(1 + ka*m);
%Large Carrier Spectrum
%%%% 6th Point (a)
figure(8)
plot(f, abs(S_LC), 'LineWidth', 0.7); 
hold on;
plot(f, abs(M_nu), 'LineWidth', 0.7);
grid on;
xlabel('frequency(hz)')
title('the spectrum of the numerical and the modulated signal')
legend('|S(f)|', '|M(f)|')
Band_Width_SLC= calc_BW(S_LC,f)
%%%% 6th Point (b)
figure(9)
plot(t, g, 'LineWidth', 0.7); 
hold on;
plot(t, s, 'LineWidth', 0.7);
grid on;
xlabel('time(second)')
title('the output of the envelope detector and am signal')
legend('g(t)', '(1 + ka*m(t)*c)')
%%% 6th Point(c)
figure(10)
plot(t, g, 'LineWidth', 0.7); 
hold on;
plot(t, m, 'LineWidth', 0.7);
grid on;
xlabel('time(second)')
title('the output of the envelope detector and the original message')
legend('g(t)', 'm(t)')
%%%% Figure 10 Spectrum 
G=fftshift(fft(g))*ts;
figure(11)
plot(f, abs(G), 'LineWidth', 0.7); 
hold on;
plot(f, abs(M_nu), 'LineWidth', 0.7);
grid on;
xlabel('Frequency (hz)')
title('the output of the envelope detector and the original message Spectrum')
legend('G(f)', 'M(f)')
%%%%%%%%%%%5 6th Point (d)%%%%%%%%%%%%%%%%%%%%
ka=1.5; % overmodulated signal 
c =cos(2*pi*fc*t);
s =(1 + ka*m).*c;                
g  =abs(hilbert(s)); %envelope detector output (demodulated signal)        
S_LC  =fftshift(fft(s))*ts;   
am =(1 + ka*m);
%Large Carrier Spectrum
%6th Point (a) Revisited
figure(12)
plot(f, abs(S_LC), 'LineWidth', 0.7); 
hold on;
plot(f, abs(M_nu), 'LineWidth', 0.7);
grid on;
xlabel('frequency(hz)')
title('the spectrum of the numerical and the modulated signal')
legend('|S(f)|', '|M(f)|')
Band_Width_SLC2= calc_BW(S_LC,f)
%%%% 6th Point (b) revisited
figure(13)
plot(t, g, 'LineWidth', 0.7); 
hold on;
plot(t, s, 'LineWidth', 0.7);
grid on;
xlabel('time(second)')
title('the output of the envelope detector and am signal')
legend('g(t)', '(1 + ka*m(t)*c)')
