function error_count5 = subblok8berptspa(seri_data3,BD_data,symbol_per_carrier,subcarrier,Phase_Set,baseband_datalength,baseband_data,GI,EbNodB,error_count);
x=seri_data3;
c1=0.90;
c3=0.40;
c5=0.17321;
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

xx=zeros(1,symbol_per_carrier*(subcarrier+GI));
for k=0:(symbol_per_carrier-1) 
    for i=1:subcarrier 
        xx(1,i+GI+(subcarrier+GI)*k)=PA(1,i+subcarrier*k);
    end 
    for i=1:GI 
        xx(1,i+(GI+subcarrier)*k)=PA(1,i+subcarrier-GI+subcarrier*k);
    end
for count=1:(length(EbNodB))
        yy=awgn(xx,EbNodB(count)+10*log10(2),'measured'); 
        for k=0:(symbol_per_carrier-1) 
            for i=1:subcarrier 
                Data_remo(1,i+subcarrier*k)=yy(1,i+GI+(GI+subcarrier)*k);
            end 
        end 
        para_data1=reshape(Data_remo,subcarrier,symbol_per_carrier);    

% ------------ FFT & DEMODULASI ----------- %
        Data_fft_after1=fft(para_data1); 
        Data_fft_after=Data_fft_after1./BD_data;  %fft pts

         for i=1:subcarrier
        for n=1:symbol_per_carrier  
            if real(Data_fft_after(i,n))>=0 
                deqpsk(i,2*n)=0; 
            else 
                deqpsk(i,2*n)=1; 
            end 
            if imag(Data_fft_after(i,n))>=0 
                deqpsk(i,2*n-1)=0; 
            else 
                deqpsk(i,2*n-1)=1; 
            end 
        end 
    end   
        Rec_data=reshape(deqpsk,1,baseband_datalength); %P/S conv
        error_count(1,count)=length(find(baseband_data~=Rec_data))+error_count(1,count); 
    end
error_count5 = error_count;  
error_count=zeros(1,length(EbNodB));
end