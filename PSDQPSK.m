clear all; 
close all; 
clc; 
bit_per_symbol=2;                   
subcarrier=1024;                  
symbol_per_carrier=100;            
% IF=4;
% N=subcarrier*symbol_per_carrier*IF;
        baseband_data=round(randint(1,symbol_per_carrier*subcarrier*2)); 
        data=reshape(baseband_data,subcarrier*2,symbol_per_carrier);
        para_data=data';
% -------- QPSK ------------

for i=1:symbol_per_carrier
   for i2=1:subcarrier
        if para_data(i,2*i2-1)==0&para_data(i,2*i2)==0 
            qpskmod(i,i2)=1+j; 
        elseif para_data(i,2*i2-1)==0&para_data(i,2*i2)==1 
            qpskmod(i,i2)=-1+j; 
        elseif para_data(i,2*i2-1)==1&para_data(i,2*i2)== 0
            qpskmod(i,i2)=-1-j; 
        else 
            qpskmod(i,i2)=1-j; 
        end 
    end
end
qpsk=qpskmod;
ofdm=ifft(qpsk,subcarrier); %IFFT ofdm, 512,symbol r+ji sebanyak 512 misal simbol 10, [512,2]

%**************************** Wiener PA ******************** 
x=ofdm;
c1 = 0.90; 
c3 = 0.40;
c5 = 0.17321

 for n = 1:length(x)
  v(n) = c1*x(n) + c3*x(n)*(abs(x(n)))^2+ c5*x(n)*(abs(x(n)))^3; 
    if n == 1 
       PA(n) =v(n);
    elseif n == 2
       PA(n)=PA(n-1)+v(n);  
    else
         PA(n)=0.2*PA(n-1)+v(n)+0.3*v(n-2); 
    end
 end

     %Predistorter Wiener
    for n = 1:length(x) 

    if n == 1 
u(n) = x(n);
    elseif n == 2
 u(n)=x(n)-0.2*u(n-1);  
    else
u(n)=-0.2*x(n-2)+x(n)-0.2*u(n-1);
    end
 PD2(n) = c1*u(n)+c3*u(n)*(abs(u(n)))^2+c5*u(n)*(abs(u(n)))^3;         
% PA Hammer
if n == 1 
       PAPD(n) =PD2(n);
    elseif n == 2
       PAPD(n)=0.2*PAPD(n-1)+PD2(n);  
    else
         PAPD(n)=0.2*PAPD(n-1)+PD2(n)+0.3*PD2(n-2);
    end
   end
% figure(9)
Fs=512;
[pxx1,f1]=pwelch(ofdm,Fs);
pxxdB1=10*log10(pxx1);
[pxx2,f2]= pwelch(PA,Fs);                  
pxxdB2=10*log10(pxx2); 
[pxx3,f3]= pwelch(PAPD,Fs);       
pxxdB3 = 10 * log10(pxx3);
%PSD 
ts=1/Fs
for k=1:Fs
    f1(k)=f1(k)/(Fs*ts);
f2(k)=f2(k)/(Fs*ts);
f3(k)=f3(k)/(Fs*ts);
end
figure(6)
plot(f1(1:Fs),[pxxdB1(Fs/2:Fs);pxxdB1(1:Fs/2-1)],'b-','linewidth',2)
hold on;
plot(f2(1:Fs),[pxxdB2(Fs/2:Fs);pxxdB2(1:Fs/2-1)],'m-','linewidth',2)
hold on;
plot(f3(1:Fs),[pxxdB3(Fs/2:Fs);pxxdB3(1:Fs/2-1)],'g-','linewidth',2)
grid on;
legend('Sinyal asli','PA','PAPD');
ylabel('dB/Hz');
xlabel('Frequency(Hz)');

% figure(7);  
% plot(f1,pxxdB1,'b','linewidth',2); 
% hold on; 
% plot(f2,pxxdB2,'r','linewidth',2); 
% hold on;
% plot(f3,pxxdB3,'g','linewidth',2); 
% hold on;
% grid on;
% legend('Sinyal asli','PA','PAPD');
% ylabel('dB/Hz');
% xlabel('Frequency(Hz)');

save 'E:\semester 8\LANJUT TA\Berkas TA\Hammer\PTS FUNGSI ENHANCED\CCDF dan BER Beda Subblok\untuk rayleigh\PSDQPSK.mat'