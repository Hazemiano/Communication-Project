clc;
close all;
clear;
ts=1/1000;  %Time step
fm=18;       %Massage Frequency
T=10/fm;      %Total simulation time
fs=0.1/ts;     % Sampling frequency
N=ceil(T/ts); %Length of Vector 
t=0:ts:(N-1)*ts;
df=fs/N;      %Frequency step
m=cos(2*pi*fm*t); %Massage Signal
m(t>9|t<0)=0;
figure (1);
plot(t,m);
xlabel('Time (s)');
ylabel('m(t)');
title('Signal m(t) = cos(2*pi*18*t)');
axis tight;
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(rem(N,2)==0)
f = - (0.5*fs) : df : (0.5*fs-df) ; 
else
    f = - (0.5*fs-0.5*df) : df : (0.5*fs-0.5*df) ; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical fft
M_nu=fftshift(fft(m))/N;


% Analytical expression %%%%% محتاجين نسال المعيد فيه %%%%
M_an = (9/2)*(sinc(9*(f-18)) + sinc(9*(f+18)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot
figure (2);
plot(f,abs(M_nu)/(max(M_nu)),"Color","b","LineWidth",1);
hold on ;
plot(f,abs(M_an)/(max(M_an)),"Color","r","LineStyle","--","LineWidth",1);
xlabel('Frequency (Hz)');
ylabel=('M(F)');
title(' Analytical and Numerical Fourier Transform');
legend('Numerical', 'Analytical');
ylim([-0.1 , 1.1]);
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
