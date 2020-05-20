function error_count1=subblok4berptspapd(PAPD,BD_data,symbol_per_carrier,subcarrier,Phase_Set,baseband_datalength,baseband_data,GI,EbNodB,error_count);
% ------------ KANAL AWGN ----------- %
xx=zeros(1,symbol_per_carrier*(subcarrier+GI)); % penambahan GI, berisi 0 wadah yang nantinya akan diisi (1,3200)
for k=0:(symbol_per_carrier-1) 
    for n=1:subcarrier  %sampai 512 kali (data subcarriernya)
        xx(1,n+GI+(subcarrier+GI)*k)=PAPD(1,n+subcarrier*k); % menjadikan Gi didepan dan dibelakang
    end 
    for n=1:GI  %i=(1 sampai 128) GI= 1/4*N=128
        xx(1,n+(GI+subcarrier)*k)=PAPD(1,n+subcarrier-GI+subcarrier*k);%penghapusan GI yg dibelakang dipindah kedepan (copyan simbol)
    end
    for count=1:(length(EbNodB))
        yy=awgn(xx,EbNodB(count)+10*log10(2),'measured');  %contoh nilai yy 0.006333875699893 + 0.019767245329202i sebanyak 640. Fungsi kanal awgn 
        for k=0:(symbol_per_carrier-1) 
            for n=1:subcarrier 
                Data_remo(1,n+subcarrier*k)=yy(1,n+GI+(GI+subcarrier)*k); %penghapusan GI
            end 
        end 
        para_data1=reshape(Data_remo,subcarrier,symbol_per_carrier);     %data yang diterima masih berupa bilangan kompleks di ubah matriksnya sebanyak 512 baris (512,1) 10 symbol 2 kolom (512,2)
% ------------ FFT & DEMODULASI ----------- %
        Data_fft_after1=fft(para_data1); 
        Data_fft_after=Data_fft_after1./BD_data;  %fft pts terdiri dari data real dan imaginer untuk pts dibagi dengan BD_data yang berisi fase terpilih yang telah dimasukkan simbol2

         for t=1:subcarrier
        for n=1:symbol_per_carrier  
            if real(Data_fft_after(t,n))>=0  %membalikkan data diambil realnya kemudian diubah menjadi 0 dan 1 lagi 
                deqpsk(t,2*n)=0; 
            else 
                deqpsk(t,2*n)=1; 
            end 
            if imag(Data_fft_after(t,n))>=0 
                deqpsk(t,2*n-1)=0; 
            else 
                deqpsk(t,2*n-1)=1; 
            end 
        end 
    end 
        Rec_data=reshape(deqpsk,1,baseband_datalength); %P/S conv
        error_count(1,count)=length(find(baseband_data~=Rec_data))+error_count(1,count); 
    end
error_count1 = error_count;  
error_count=zeros(1,length(EbNodB));
end
