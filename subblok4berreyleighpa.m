function error_count4=subblok4berreyleighpa(PA,BD_data,symbol_per_carrier,subcarrier,Phase_Set,baseband_datalength,baseband_data,GI,EbNodB,error_count);

% ------------ KANAL reyleigh ----------- %
xx=zeros(1,symbol_per_carrier*(subcarrier+GI)); %data ditampilkan seri, dgn matrik[ 1, 1*((512+(1/4*512))]= [1, 640]

   for k=0:(symbol_per_carrier-1) 
    for i=1:subcarrier  %sampai 512 kali (data subcarriernya)
        xx(1,i+GI+(subcarrier+GI)*k)=PA(1,i+subcarrier*k); % misal 1 data nilainya = 0.019531250000000 - 0.019531250000000i
    end 
    for i=1:GI  %i=(1 sampai 128) GI= 1/4*N=128
        xx(1,i+(GI+subcarrier)*k)=PA(1,i+subcarrier-GI+subcarrier*k); %imaginer seri
    end

        for count=1:(length(EbNodB))    
           snr= EbNodB(count) + 10*log10(2);
             sgma=2/sqrt(snr*2); %variance
             u=rand;% uniform random variable in (0,1)
             z=sgma*(sqrt(2*log(1/(1-u)))); % rayleigh distributed random variable
          
          % rayleigh function 
           uu=2*pi*(xx);
           n=z+(sqrt(sgma/1))*rand(1,1);
           xy= n.*exp(uu);
           yy=awgn(xy,snr,'measured');
for k=0:(symbol_per_carrier-1) 
                for u=1:subcarrier
                    Data_remo(1,u+subcarrier*k)=yy(1,u+GI+(GI+subcarrier)*k);
                end 
            end 
            para_data1=reshape(Data_remo,subcarrier,symbol_per_carrier);
% ------------ FFT & DEMODULASI ----------- %
        Data_fft_after1=fft(para_data1); 
        Data_fft_after=Data_fft_after1./BD_data;  %fft pts

         for u=1:subcarrier
        for n=1:symbol_per_carrier  
            if real(Data_fft_after(u,n))>=0 
                deqpsk(u,2*n)=0; 
            else 
                deqpsk(u,2*n)=1; 
            end 
            if imag(Data_fft_after(u,n))>=0 
                deqpsk(u,2*n-1)=0; 
            else 
                deqpsk(u,2*n-1)=1; 
            end 
        end 
    end 
        Rec_data=reshape(deqpsk,1,baseband_datalength); %P/S conv
        error_count(1,count)=length(find(baseband_data~=Rec_data))+error_count(1,count); 
    end
error_count4 = error_count;  
error_count=zeros(1,length(EbNodB));
end
end