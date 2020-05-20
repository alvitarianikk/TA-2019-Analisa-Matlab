function [PAPD8]=PAPD8(seri_data3);

x = seri_data3;

%koefisien 

c1=0.90;
c3=0.40;
c5=0.17321;
% ======================= LINEARISASI =====================%
% ---------------- PD Hamm ---------------- %
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
PA=PA;
     %Predistorter Wiener
    for n = 1:length(x) 

    if n == 1 
u(n) = x(n);
    elseif n == 2
 u(n)=x(n)-0.2*u(n-1);  
    else
u(n)=-0.3*x(n-2)+x(n)-0.2*u(n-1);
    end
 PD(n) = c1*u(n)+c3*u(n)*(abs(u(n)))^2+c5*u(n)*(abs(u(n)))^3;      
   PD=PD;
if n == 1 
       PAPD(n) =PD(n);
    elseif n == 2
       PAPD(n)=0.2*PAPD(n-1)+PD(n);  
    else
         PAPD(n)=0.2*PAPD(n-1)+PD(n)+0.3*PD(n-2);
end
  PAPD8=PAPD;
   end 
