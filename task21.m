clc
clear all
ZL=50;
Z0=100;

N=(1:5); 
Z = zeros(1,length(N));
Z = [Z0 Z ZL];
hold all; 
theta=(0:0.01:pi); 
df=3-4/pi*theta; 

for n = 1:length(N) 
    cn=nchoosek(length(N),n-1);
    Z(n+1) = exp(log(Z(n))+2^(-1*length(N))*cn*log(ZL/Z0));
end

 for n=N 
     A=1/2^n*(ZL-Z0)/(ZL+Z0);  
     g2=2^n*abs(A)*abs((cos(theta)).^n);
     plot(df,g2); 
 end
 legend('N=1','N=2','N=3','N=4','N=5')















