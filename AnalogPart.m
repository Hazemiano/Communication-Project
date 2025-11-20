clc;
close all;
clear;
ts=1/600;  %Time step
fm=18;       %Massage Frequency
T=10/fm;      %Total simulation time
fs=1/ts;     % Sampling frequency
N=ceil(T/ts); %Length of Vector 
t=0:ts:((N-1)*ts);
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
plot(f,abs(M_nu),"Color","b","LineWidth",1);
hold on ;
plot(f,abs(M_an),"Color","r","LineStyle","--","LineWidth",1);
xlabel('Frequency (Hz)');
ylabel=('|M(F)|');
title(' Analytical and Numerical Fourier Transform');
legend('Numerical', 'Analytical');
ylim([-0.1 , 1.1]);
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%get Band_Width
Max_Value = max(abs(M_nu));
index = find(f == 0);
b = 0;

for C_index = index : length(f)


    if (abs(M_nu(C_index)) == Max_Value)
        b = 1;
    end

    if(M_nu(C_index) < (0.01 *  Max_Value) && b == 1)
    index = C_index;
    break
    end
end


Band_Width = f(index);

disp(['Band Width = ' num2str(Band_Width)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SSB Modulation (S & s)
fc = 50;
Ac = 10;
c = Ac * cos(2*pi*fc*t);
s = m.* c;
S = fftshift(fft(s))/N; %%%%% S is after DSB-SC then it goes to BPF to take wanted side

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1st method by BPF
BPF = zeros(size(f));
BPF(f>(fc-Band_Width) & f<(fc))=1;%% +ve frequency
BPF(f<-(fc-Band_Width) & f>-(fc))=1;%% -ve frequency
S_LSB1 = BPF.*S;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2nd method by LBF 
LBF = abs(f) <fc; %%%%LBF
S_LSB2 = LBF.* S;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure (3);
plot(f,abs(S_LSB2),"Color","b" , "LineWidth",0.6);
hold on ;

xlabel('Frequency (Hz)');
ylabel=('|S(F)|_LBS2');
title('Modulated Message 2St Method');


plot(f,abs(S_LSB1),"Color","r");
hold on ;
grid on;
xlabel('Frequency (Hz)');
ylabel=('|S(F)|_LBS1');
title('Modulated Message 1St Method');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% S1_Band_Width (BW for S_LSB1)
Max_Value = max(abs(S_LSB1));
index = find(f == 0);
b = 0;

for C_index = index : length(f)


    if (abs(S_LSB1(C_index)) == Max_Value)
        b = 1;
    end

    if(S_LSB1(C_index) < (0.01 *  Max_Value) && b == 1)
    index = C_index;
    break
    end
end


S1_Band_Width = f(index);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  S2_Band_Width (BW for S_LSB2)
Max_Value = max(abs(S_LSB2));
index = find(f == 0);
b = 0;

for C_index = index : length(f)


    if (abs(S_LSB2(C_index)) == Max_Value)
        b = 1;
    end

    if(S_LSB2(C_index) < (0.01 *  Max_Value) && b == 1)
    index = C_index;
    break
    end
end


S2_Band_Width = f(index);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% print BW of 2 mehtods
disp(['1st method BW_SSB = ' num2str(S1_Band_Width)]);
disp(['2nd method BW_SSB = ' num2str(S2_Band_Width)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End point 5A 
r=ifft(fftshift(S_LSB1))*N;
v=r .* c;
V=fftshift(fft(v))/N;
H=abs(f)<18;

mo=real(ifft(fftshift(H.* V))*N);
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





