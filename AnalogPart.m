clc;
close all;
clear;
fc = 50;
Ac = 1;
ts=0.01/fc;  %Time step
fm=18;       %Message Frequency
T=10/fm;      %Total simulation time
fs=1/ts;     % Sampling frequency
N=ceil(T/ts); %Length of Vector 
t=0:ts:(N-1)*ts;
df=fs/N;      %Frequency step
m=cos(2*pi*fm*t); %Message Signal
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical fft
M_nu=fftshift(fft(m))/N;


% Analytical expression 
M_an = 4.5*sinc(9*(f-18)).*exp(-i*pi*9*(f-18)) + 4.5*sinc(9*(f+18)).*exp(-i*pi*9*(f+18));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot
figure (2);
plot(f,abs(M_nu),"Color","b","LineWidth",1);
hold on ;
plot(f,abs(M_an),"Color","r","LineStyle","--","LineWidth",1);
xlabel('Frequency (Hz)');
ylabel ('M(F)');
title(' Analytical and Numerical Fourier Transform');
legend('Numerical', 'Analytical');
axis tight;
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% get Band_Width
Max_Value = max(abs(M_nu));
Therhold_Value = 0.01 * Max_Value;
Above_Therhold = abs (M_nu) >= abs(Therhold_Value);

%figure (88)
%plot (f,Above_Therhold)

Freq_Index = find(Above_Therhold);

if ~isempty(Freq_Index)
    % These define the bandwidth higher edge
    f_max = f(Freq_Index(end));
       
end
    
Band_Width = f_max;

disp(['Band Width = ' num2str(Band_Width)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Point 5
c = Ac * cos(2*pi*fc*t);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1st method by BPF
s=m.*c;

S=fftshift(fft(s))/N; % Modulated Signal in freq domain
figure(4) 
plot(f,abs(S));
hold on ;
xlabel('Frequency (Hz)')
ylabel ('S(F)');
title('Modulated Signal 5A' )
grid on;
axis tight;

BPF = zeros(size(f));
BPF(f>(fc-Band_Width) & f<=(fc))=1;    % +ve frequency
BPF(f<-(fc-Band_Width) & f>= -(fc))=1; % -ve frequency
S_LSB1 = BPF.*S;  % 1st LSB Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2nd method by Hilbert transform of m(t)
m_hilbert = imag(hilbert(m));

% LSB modulation: s(t) = m(t)*cos(wc*t) + m_hat(t)*sin(wc*t)
s_LSB2 = 0.5 * m .*  cos(2*pi*fc*t) + 0.5 *m_hilbert.* sin(2*pi*fc*t);

% LSB in freq doman (method 2)
S_LSB2 = fftshift(fft(s_LSB2)/N); %2nd LSB Signal in freq domain 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot of S_LSB1 & S_LSB2
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% S1_Band_Width (BW for S_LSB1)


Max_Value = max(abs(S_LSB1));
Therhold_Value = 0.01 * Max_Value;
Above_Therhold = abs (M_nu) > abs(Therhold_Value);
Freq_Index = find(Above_Therhold);

if ~isempty(Freq_Index)
    % These define the bandwidth edges
    f_min = f(Freq_Index(1));
    f_max = f(Freq_Index(end));
       
end
    

S1_Band_Width = f_max - f_min;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  S2_Band_Width (BW for S_LSB2)
Max_Value = max(abs(S_LSB2));
Therhold_Value = 0.01 * Max_Value;
Above_Therhold = abs (M_nu) >= abs(Therhold_Value);
Freq_Index = find(Above_Therhold);

if ~isempty(Freq_Index)
    % These define the bandwidth edges
    f_min = f(Freq_Index(1));
    f_max = f(Freq_Index(end));
       
end
    
S2_Band_Width = f_max - f_min;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% print BW of 2 mehtods
disp(['1st method BW_SSB = ' num2str(S1_Band_Width)]);
disp(['2nd method BW_SSB = ' num2str(S2_Band_Width)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End point 5A 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% demodulated signal Using Coherent Detector


s_LSB2 = ifft(fftshift(S_LSB2))*N; % LSB in time domain 
v = s_LSB2.*c;

V=fftshift(fft(v))/N; % Demodulated Signal before LBF

H=abs(f)<22; %LBF 

mo=real(ifft(fftshift(H.* V))*N); % Demodilated Signal

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


