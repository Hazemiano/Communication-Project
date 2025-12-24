clc;
%pkg load signal
close all;
clear;
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
[maxi, index] = max(M_nu);
threshold = 0.01 * maxi;
threshold_crossing = M_nu - threshold;
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
Band_Width=21.6;
BPF = zeros(size(f));
BPF(f>(fc-Band_Width) & f<(fc))=1; %% +ve frequency
BPF(f<-(fc-Band_Width) & f>-(fc))=1;%% -ve frequency
S =  S.*BPF;
S_LSB1=S; 
%%%%%%%%%%%$$$$$$$$ 2nd method
LPF = abs(f) <fc; %%%%LPF
S_LSB2 = LPF.* S;
S= S.*LPF;
S_LSB2 =S;

figure (5);
plot(f,abs(S_LSB1),"Color","r");
hold on ;
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



