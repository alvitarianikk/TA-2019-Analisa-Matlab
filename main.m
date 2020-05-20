close all; 
clc; 
format long;
clear all; 
%------------------------PROGRAM UTAMA------------------------------------%

%=============================Komponen PTS================================%
bit_per_symbol=2;               %jumlah bit qpsk
subcarrier=512;               %jumlah subcarrier (N)
symbol_per_carrier=1000;          %jumlah simbol OFDM                      
PrefixRatio=1/4;                %secara umum, rentang CP dalam sistem OFDM berkisar dari 1/4 hingga 1/32 dari periode simbol.
GI=PrefixRatio*subcarrier;    %adalah rasio waktu Cyclic Prefix "CP" ke inverse FFT "T(IFFT). GI digunakan untuk menghilangkan inter-simbol dan inter-carrier interferensi.
                                %Salinan GI terakhir T (GI) dari periode simbol yang berguna "T (IFFT)",
                                %disebut Cyclic Prefix "CP", digunakan untuk mengumpulkan banyak, sambil mempertahankan ortogonalitas dari subcarrier.
                                %Setiap simbol ditransmisikan untuk waktu yang sedikit lebih lama, waktu simbol yang diperpanjang T (s), daripada waktu simbol yang aktif (atau berguna) T (IFFT).
EbNodB=1:1:30;       
Eb=sqrt(2)/2;   
Phase_Set = [1 -1]; 

%==================Pembangkitan bit & modulasi============================%
BD_data=zeros(subcarrier,symbol_per_carrier);
error_count=zeros(1,length(EbNodB));
baseband_datalength=subcarrier*symbol_per_carrier*bit_per_symbol;

for z=1:3 %untuk kondisi ofdm, pts4 dan pts8 enhanced 
   Rec_data=zeros(1,baseband_datalength); 
   baseband_data=randint(1,baseband_datalength,2); %pembangkitan biner 01 sebanyak baseband data length
   para_data=reshape(baseband_data,subcarrier,symbol_per_carrier*bit_per_symbol); %S/P converter
 % ------------------------ QPSK --------------------------------- %
data_length=baseband_datalength/subcarrier;
for i=1:subcarrier
   for i2=1:data_length/bit_per_symbol
        if para_data(i,2*i2-1)==0&para_data(i,2*i2)==0 
            qpskmod(i,i2)=1+j; 
        elseif para_data(i,2*i2-1)==0&para_data(i,2*i2)==1 
            qpskmod(i,i2)=-1+j; 
        elseif para_data(i,2*i2-1)==1&para_data(i,2*i2)==1 
            qpskmod(i,i2)=-1-j; 
        else 
            qpskmod(i,i2)=1-j; 
        end 
    end
end
          qpsk=qpskmod;%128*100
% ---------------OFDM ---------- %   
        after_qpsk=ifft(qpsk,subcarrier); %IFFT ofdm, 512,symbol r+ji sebanyak 512 misal simbol 10, [512,2]
        ofdm = reshape(after_qpsk,[],subcarrier);
        Signal_Power =abs(ofdm.^2);
        Peak_Power = max(Signal_Power,[],2);
        Mean_Power = mean(Signal_Power,2);
        PAPR_Original = 10*log10(Peak_Power./Mean_Power); %papr ofdm
        PAPR_Original1= min(PAPR_Original); 
       
       if(z == 1)
            seri_data1=reshape(after_qpsk,1,symbol_per_carrier*subcarrier); %PS converter
       elseif(z == 2)
            [cdf2, PAPR2,seri_data2,BD_data]=subblok4ccdf(qpsk,BD_data,symbol_per_carrier,subcarrier,Phase_Set);
            [out1,out2,out3]=PAPD(seri_data2);
        x=seri_data2;
        PA=out1;
        PD=out2;
        PAPD=out3;
            error_count1 = subblok4berptspapd(PAPD,BD_data,symbol_per_carrier,subcarrier,Phase_Set,baseband_datalength,baseband_data,GI,EbNodB,error_count);
            error_count2 = subblok4berptspa(PA,BD_data,symbol_per_carrier,subcarrier,Phase_Set,baseband_datalength,baseband_data,GI,EbNodB,error_count);
            error_count3 = subblok4berreyleigh(PAPD,BD_data,symbol_per_carrier,subcarrier,Phase_Set,baseband_datalength,baseband_data,GI,EbNodB,error_count);
            error_count4 = subblok4berreyleighpa(PA,BD_data,symbol_per_carrier,subcarrier,Phase_Set,baseband_datalength,baseband_data,GI,EbNodB,error_count);
       else
           [cdf3, PAPR3,seri_data3,BD_data]=subblok8ccdf(qpsk,BD_data,symbol_per_carrier,subcarrier,Phase_Set);
           error_count5 = subblok8berptspa(seri_data3,BD_data,symbol_per_carrier,subcarrier,Phase_Set,baseband_datalength,baseband_data,GI,EbNodB,error_count);
       error_count6 = subblok8berptspapd(seri_data3,BD_data,symbol_per_carrier,subcarrier,Phase_Set,baseband_datalength,baseband_data,GI,EbNodB,error_count);


       end 
end
% %==================Plot Hasil============================%
[cdf1, PAPR1] = ecdf(PAPR_Original);
%PLOT CCDF dalam dB
figure(1)
semilogy(PAPR1,1-cdf1,'*-g',PAPR2,1-cdf2,'b',PAPR3,1-cdf3,'m');
grid on
hold on
legend('OFDM Asli','PTS Subblok 4','PTS Subblok 8')
title ('Kurva CCDF');
xlabel('PAPR0 [dB]'); 
ylabel('CCDF (Pr[PAPR>PAPR0])');

% %PLOT AM/AM
% figure (2)
% plot(abs(x),abs(PA),'b.');
% hold on;
% grid on;
% title('PA')
% xlabel('vin (volt)');
% ylabel('vout (volt)');
% figure (3)
% plot(abs(x),abs(PD),'r.');
% hold on;
% grid on; 
% title('PD')
% xlabel('vin (volt)');
% ylabel('vout (volt)');
% figure(4)
% plot(abs(x),abs(PAPD),'g.');
% hold on;
% grid on; 
% title('PAPD')
% xlabel('vin (volt)');
% ylabel('vout (volt)');
% % 
% figure (5)
% plot(abs(x),abs(PA),'b.');
% hold on;
% grid on; 
% plot(abs(x),abs(PD),'r.');
% hold on;
% grid on; 
% plot(abs(x),abs(PAPD),'g.');
% hold on;
% grid on; 
% title('Kurva AM/AM Power Amplifier');
% legend('PA model', 'PD model', 'linearisasi');
% xlabel('vin (volt)');
% ylabel('vout (volt)');

%Hitung nilai BER terhadap SNR
error_count11=error_count1./(baseband_datalength)
error_count12=error_count2./(baseband_datalength)
error_count13=error_count3./(baseband_datalength)
error_count14=error_count4./(baseband_datalength)
error_count15=error_count5./(baseband_datalength)
error_count16=error_count6./(baseband_datalength)

% %PLOT BER 
figure(6)
semilogy(EbNodB,error_count11,'k+-','linewidth',2); 
hold on; 
semilogy(EbNodB,error_count12,'b+-','linewidth',2); 
hold on; 
semilogy(EbNodB,error_count13,'r+-','linewidth',2); 
hold on; 
semilogy(EbNodB,error_count14,'g+-','linewidth',2); 
hold on;
legend('PTS PAPD  Subblok 4','PTS PA Subblok 4','PTS PAPD Reyleigh V=4','PTS PA Reyleigh V=4')
title ('Kurva BER');
xlabel('SNRdb(dB)');
ylabel('BER');
axis([0 30 10^(-7) 10^(0)]); 
hold on 
grid on;

figure(7)
semilogy(EbNodB,error_count11,'k+-','linewidth',2); 
hold on; 
semilogy(EbNodB,error_count12,'b+-','linewidth',2); 
hold on; 
semilogy(EbNodB,error_count15,'m+-','linewidth',2); 
hold on;
semilogy(EbNodB,error_count16,'y+-','linewidth',2); 
hold on;
legend('PTS PAPD  Subblok 4','PTS PA Subblok 4','PTS PA V=8','PTS PAPD V=8')
title ('Kurva BER');
xlabel('SNRdb(dB)');
ylabel('BER');
axis([0 30 10^(-7) 10^(0)]); 
hold on 
grid on;
save 'E:\semester 8\LANJUT TA\Berkas TA\Hammer\PTS FUNGSI ENHANCED\CCDF dan BER Beda Subblok\untuk rayleigh/main.mat'
% save 'E:\semester 8\LANJUT TA\Berkas TA\Hammer\PTS FUNGSI ENHANCED\CCDF dan BER Beda Subblok\untuk rayleigh/mainam.mat'