function [cdf2, PAPR2,seri_data2,BD_data]=subblok4ccdf(qpsk,BD_data,symbol_per_carrier,subcarrier,Phase_Set);

Vt = 4;   
V= 2*Vt; %Perbedaan PTS enhanced dgn biasa, terdapat peningkatan jumlah subblok, namun kombinasi tetap. N/2V
Choose = [1 1 1 1 1 1 1 1 ; 1 1 1 2 2 1 1 1 ; 1 1 2 1 1 2 1 1; 1 2 1 1 1 1 2 1; 2 1 1 1 1 1 1 2;... 
          1 1 2 2 2 2 1 1 ; 1 2 1 2 2 1 2 1 ; 1 2 2 1 1 2 2 1; 2 2 1 1 1 1 2 2; 2 1 2 1 1 2 1 2 ;...
          2 1 1 2 2 1 1 2 ; 2 2 2 1 1 2 2 2 ; 2 2 1 2 2 1 2 2; 2 1 2 2 2 2 1 2; 1 2 2 2 2 2 2 1;...
          2 2 2 2 2 2 2 1];    
Choose_Len = 16;   
                           % size 16x(4x2) % kombinasi faktor fase dengan jumlah subblok 4 yg di enhanced

qpsk=qpsk.'; %P/S
BD_data=BD_data'; 
for ii=1:symbol_per_carrier     
    A = zeros(V,subcarrier); 
  for k2=1:V  
        A(k2,k2:V:subcarrier)=qpsk(ii,k2:V:subcarrier);  %partisi subblok dengan partisi interleaved
  end   
    a = ifft(A,[],2);

% ------------ PTS ----------- 
  min_value = 10;  %treshold terserah
  for n=1:Choose_Len 
        temp_phase = Phase_Set(Choose(n,:)).'; %mengubah kombinasi yg tadinya misal 1212 menjadi 1 -1 1 -1 ditransposes dari bentuk seri ke paralel
        temp_max = max(abs(sum(a.*repmat(temp_phase,1,subcarrier))));  %perkalian dengan faktor fasa repmat untuk mengcopy faktor fasa ke seluruh sublok. Lalu diiterasi dgn mengganti dengan faktor fasa lain untuk diproses selanjutnya 

     if temp_max<min_value            % apabila puncak maksimum < threshold, maka puncak maksimum lah terkecil %sinyal kandidat dengan nilai terkecil
            min_value = temp_max;           %16 kombinasi faktor fase diproses satu per satu dicari nilai perkalian (temp_max) yang paling kecil 12121111 11112121
            Best_n = n;                     %urutan kebrapa yang paling kecil nilainya? misal yg paling bagus adalah kombinasi ke 5. Lalu looping hingga kombinasi paling akhir. Misal ditemukan kombinasi 6
     end 
  end 
    qpsk(ii,:) = sum(a.*repmat(Phase_Set(Choose(Best_n,:)).',1,subcarrier));  %fasa terbaik tadi dicopy (fungsi repmat) sebanyak (simbol,subcarrier) dikalikan sinyal a
    use_phase=Phase_Set(Choose(Best_n,:)); %mengubah kombinasi lagi misal yg tadinya paling bagus adalah 1211 menjadi 1 -1 1 1
    Signal_Power4 =abs(qpsk.^2);
    Peak_Power4 = max(Signal_Power4,[],2);
    Mean_Power4 = mean(Signal_Power4,2);
    PAPR_PTS4= 10*log10(Peak_Power4./Mean_Power4);
  for k2=1:V 
        no=repmat(use_phase(k2),1,subcarrier);  % hasil dari use phase diinisialisasi sebagai k2[1,512]berisi 0 sebanyak subbloknya (8)
        BD_data(ii,k2:V:subcarrier)=no(k2:V:subcarrier); %side information. Mengisi dengan fase terpilih (use_phase) misal 1 -1 1 1 secara bertahap sesuai urutan partisinya (interleaved) diiterasi sebanyak simbol
  end 
end

qpsk=qpsk.';%P/S
BD_data=BD_data';%num_subcarr*100 
seri_data=reshape(qpsk,1,symbol_per_carrier*subcarrier); 
seri_data2 = seri_data;
[cdf2, PAPR2] = ecdf(PAPR_PTS4);


