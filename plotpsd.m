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

 c=(abs(PA));
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
N=subcarrier;
pxx1=pwelch(x,[],[],N);
pxxdB1=10*log10(pxx1);
pxx2= pwelch(PA,[],[],N);                  
pxxdB2=10*log10(pxx2); 
pxx3= pwelch(PAPD,[],[],N);       
pxxdB3 = 10 * log10(pxx3);
% PSD 
ts=1/8000;
for k=1:N
    fr(k)=(k-1-(N/2))/(N*ts);
end
figure(6)
plot(fr(1:N),[pxxdB1(N/2:N);pxxdB1(1:N/2-1)],'b-')
hold on;
plot(fr(1:N),[pxxdB2(N/2:N);pxxdB2(1:N/2-1)],'m-')
hold on;
plot(fr(1:N),[pxxdB3(N/2:N);pxxdB3(1:N/2-1)],'g-')
grid on;
legend('SinyDal asli','PA','PAPD');
ylabel('dB/Hz');
xlabel('Frequency(Hz)');
save('plotpsd.mat','pxxdB1','pxxdB2','pxxdB3','x','subcarrier')