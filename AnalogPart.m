clc;
close all;
clear;
fc = 50;
ts=0.1/fc;  %Time step
fm=18;       %Massage Frequency

T=10/fm;      %Total simulation time
fs=1/ts;     % Sampling frequency
N=ceil(T/ts); %Length of Vector 
t=0:ts:(N-1)*ts;
df=fs/N;      %Frequency step
m=cos(2*pi*fm*t); %Massage Signal
m(t>=9|t<=0)=0;
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
M_nu=fftshift(fft(m))*ts;


% Analytical expression 
M_an = 4.5*sinc(9*(f-18)).*exp(-i*pi*9*(f-18)) + 4.5*sinc(9*(f+18)).*exp(-i*pi*9*(f+18));

%M_an = (9/2)*(sinc(9*(f-18)) + sinc(9*(f+18)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot
figure (2);
plot(f,abs(M_nu),"Color","b","LineWidth",1);
hold on ;
plot(f,abs(M_an),"Color","r","LineStyle","--","LineWidth",1);
xlabel('Frequency (Hz)');
ylabel=('M(F)');
title(' Analytical and Numerical Fourier Transform');
legend('Numerical', 'Analytical');

grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%get Band_Width
% Assume M_nu = spectrum magnitude, f = frequency axis
[maxi, index] = max(M_nu);        
threshold = 0.01 * maxi;         

stopindex_right = index;
for i = index:length(f)
    if M_nu(i) < threshold
        stopindex_right = i;
        break
    end
end

stopindex_left = index;
for i = index:-1:1
    if M_nu(i) < threshold
        stopindex_left = i;
        break
    end
end

Band_Width = f(stopindex_right) - f(stopindex_left);

disp(['Estimated Bandwidth = ' num2str(Band_Width) ' Hz']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Ac = 10;
c = Ac * cos(2*pi*fc*t);
s = m.* c;
S = fftshift(fft(s))*ts; %%%%% S is after DSB-SC then it goes to BPF to take wanted side

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1st method by BPF
BPF = zeros(size(f));
BPF(f>(fc-Band_Width) & f<(fc))=1;%% +ve frequency
BPF(f<-(fc-Band_Width) & f>-(fc))=1;%% -ve frequency
S =  S.*BPF;
S_LSB1=S;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2nd method by LBF 
LBF = abs(f) <fc; %%%%LBF
S= S.*LBF;
S_LSB2 =S;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure (3);
plot(f,abs(S_LSB2),"Color","b" , "LineWidth",0.6);
hold on ;

xlabel('Frequency (Hz)');
ylabel=('|S(F)|_LBS2');
title('Modulated Message 2St Method');


plot(f,abs(S_LSB1),"Color","r","LineWidth",0.6);
hold on ;
grid on;
xlabel('Frequency (Hz)');
ylabel=('|S(F)|_LBS1');
title('Modulated Message 1St Method');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% S1_Band_Width (BW for S_LSB1)
[maxi, index] = max(abs(S_LSB1));
threshold = 0.01 * maxi;

stopindex_left = index;
for i = index:-1:1
    if abs(S_LSB1(i)) < threshold
        stopindex_left = i;
        break
    end
end

S1_Band_Width = f(index) - f(stopindex_left);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  S2_Band_Width (BW for S_LSB2)
[maxi, index] = max(abs(S_LSB2));
threshold = 0.01 * maxi;

stopindex_left = index;
for i = index:-1:1
    if abs(S_LSB2(i)) < threshold
        stopindex_left = i;
        break
    end
end

S2_Band_Width = f(index) - f(stopindex_left);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% print BW of 2 mehtods
disp(['1st method BW_SSB = ' num2str(S1_Band_Width)]);
disp(['2nd method BW_SSB = ' num2str(S2_Band_Width)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End point 5A 
r=ifft(fftshift(S_LSB1))/ts;
v=r .* c;
V=fftshift(fft(v))*ts;
H=abs(f)<18;

mo=real(ifft(fftshift(H.* V))/ts);
figure(4)
plot(t,mo);
hold on;
plot(t,m);
legend('demodulated massage signal','original massage')
xlabel('Time (s)');
%ylabel('m(t)');
title(' original massage & demodulated massage signal')
axis tight;
grid on;





