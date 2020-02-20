clc
clear all
%In the first graph we assume Z0=50 and Zl=100, the impendace Z is defined
%by the formula: Z=Z0*exp(a*z) where z is the length of the line we're
%calculating and a is a defined by (1/L)*log(Zl/Z0). Then we will plot the
%vector Z according to the length of the line.
L = 10; %assume length of L = 10
beta = 1;
z = 0.01:0.01:L;
Z0 = 50;
Zl = 100;
a = (1/L)*log(Zl/Z0);
Z = Z0*exp(a*z);
figure;
plot(z,Z);
ylim([40 110]);
title('Impendace to length relation');
xlabel('Length of line (z)');
ylabel('Impendace Z(z)');

%In the next graph we will have the Z0 and Zl reversed, now Z0=100 and
%Zl=50, we will use the formula abs(0.5*log(Zl/Z0)*(sin(BL(i))/(BL(i)))) to
%find Gamma for every number of sections, for example, if we want to find
%Gamma0 we will say Bl is equal to 0 and we will see the reflection 
%coefficient is actually 0.346, so basically we made a vector of N=1000 and
%found the reflection coefficient for each section and plot it on the
%graph.
Z0 = 100;
Zl = 50;
BL = linspace(0,10*pi,1000); %to calculate 1000 sections with a length of 10pi/1000
%BL = 0:0.01:10*pi; to calculate for x section with a length of 0.01
Rtheta = zeros(1,length(BL));
for i = 1:length(BL)
   Rtheta(i) = abs(0.5*log(Zl/Z0)*(sin(BL(i))/(BL(i)))); 
end
figure;
plot(BL,Rtheta);
xlabel('BL - Length of the line');
ylabel('Amplitude of reflection coefficient');