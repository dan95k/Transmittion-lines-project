%This code is used only for up to 4 sections transformer (as requested), 
%it is possible to use it to any order of chebyshevs polinomial
%just need to add its higher order polinomials to the code, but the idea is
%clear.
clc
clear all
Zl = 100;
Z0 = 50;
N = 1:3;
Rm = 0.05;
theta=(0:0.01:pi); 
secThetaM = cosh((1/length(N))*acosh((log(Zl/Z0)/(2*Rm))));
thetaM = asec(secThetaM)*180/pi;
A = ((Zl-Z0)/(Zl+Z0))*1./(chebyshevT(length(N),sec(thetaM))); %equals Gamma0
%I tried to map every coefficient to a new cell, i've used the formulas for
%each chebyshev polinomial, and for each coefficient used a new cell, for
%example if i need an order 4 chebyshev polinomial, so according to the
%formula i have cos(4x)*sec(thetaM)^4, so i placed sec(thetaM)^4 into the
%first place of the vector, now i know the coefficient of cos(4x). Next i
%know that the coefficient of cos(2x) is 4*sec(thetaM)^4-4sec(thetaM)^2 and
%placed it into the second place of the vector and etc. then i would be
%able to divide the coefficients of both sides and find the right Gammas.
%ofcourse this method require to determine each polinomial so i made it
%upto N=4 but the idea is the same for higher order polinomials.
if length(N) == 1 % Tn(cos(x)*sec(y)) = cos(x)*sec(y)
    Rside(1) = sec(thetaM)*cos(x);
    Rside = A*Rside;
elseif length(N) == 2 %Tn(cos(x)*sec(y)) = 2(cos(x)*sec(y))^2-1
    Rside(1) = (sec(thetaM))^2; Rside(2) = ((sec(thetaM))^2)+1;
    Rside = A*Rside;
elseif length(N) == 3 %Tn(cos(x)*sec(y)) = 4(cos(x)*sec(y))^3-3(cos(x)*sec(y))
    Rside(1) = (sec(thetaM))^3; Rside(2) = 3*((sec(thetaM))^3)-3*sec(thetaM);
    Rside = A*Rside;
elseif length(N) == 4 %Tn(cos(x)*sec(y)) = 8(cos(x)*sec(y))^4-8(cos(x)*sec(y))^2+1
    Rside(1) = (sec(thetaM)^4); Rside(2) = 4*(sec(thetaM)^4)-4*(sec(thetaM)^2); Rside(3) = 3*(sec(thetaM)^4)-4*(sec(thetaM)^2)+1;
    Rside = A*Rside;
end
Gammas = zeros(1,length(N)); %placeholder matrix for reflection coefficients
%We are running this loop 2 times because we know Gamma3 == Gamma0 and
%Gamma1 == Gamma2, so no need to run it N times.
for i = 1:2 
    Gammas(i) = Rside(i)/2;
end

if length(N) >= 3
   Gammas(3) = Gammas(2);
   if length(N) >= 4
    Gammas(4) = Gammas(1);
   end
end

Zmat = [Z0 ones(1,length(N)) Zl];
%Now we will calculate the characteristic impedances using the formula
%ln(Z1) = ln(Z0)+2Gamma0
j=1;
for i = 1:length(N)
    Zmat(i+1) = exp(log(Zmat(i))+2*Gammas(j));
    j = j + 1;
end
%f/f0 is nonmalized freq, its corresponding theta formula is
%2*theta/pi,
f=2*theta/pi;
figure()
hold on;
Gamma1 = abs(A*exp(-1j)*chebyshevT(1,cos(theta)*sec(thetaM)));
plot(f,Gamma1,'Linewidth',2);
legend('N=1');
if length(N) >=2
    Gamma2 = abs(A*exp(-1j*2)*chebyshevT(2,cos(theta)*sec(thetaM)));
    plot(f,Gamma2,'Linewidth',2);
    legend('N=1','N=2');
    if length(N) >= 3
        Gamma3 = abs(A*exp(-1j*3)*chebyshevT(3,cos(theta)*sec(thetaM)));
        plot(f,Gamma3,'Linewidth',2);
        legend('N=1','N=2','N=3');
        if length(N) == 4
            Gamma4 = abs(A*exp(-1j*4)*chebyshevT(4,cos(theta)*sec(thetaM)));
            plot(f,Gamma4,'Linewidth',2);
            legend('N=1','N=2','N=3','N=4');
        end
    end
end
title('|Gamma(\theta)| to frequency relation');
xlabel('Normalized Frequency');
ylabel('|Gamma(\theta)|');
xlim([1/3 5/3]);