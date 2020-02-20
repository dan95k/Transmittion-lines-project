% In this code we perform an analysis of a light passing through an array
% of layers in a window. We shoot a ray of light into the window and we
% want it to pass th window without it breaking. the ranges of the incident
% angle and wavelength as follows: from -45 degrees to 45 degrees and for
% the wavelength are the whole visible spectrum. in order to do that we use
% an anlogy to a transmittion line when we want to match the input
% impadence to the load (input impadence is the space before the window and
% the load is the space after the window) when every layer we add is
% another transmittion line. Epsilons of every metter is given.

% The idea is as follows, we calculate each layer zeta using
% the epsilons given using the forumla sqrt(?0*?r/?0?r).
% Later we find each layer width as the width changes for every wavelngth
% and angle of incidance. Beta is calculated with the forumla 2*pi*f(i)*n(k)/c
% when f is the frequency and n is the refraction index of the layer and c
% is the speed of light.
% As in transmittion lines we have an characteristic impendace so we have
% here, we calculate it here (zeta of the metter)/cos(the angle in the layer)
% after we have the characteristic impendaces beta and the width of every
% layer we can procced to calculate the reflection coefficients.
% We use the formula Z0*( ZL+jZ0*tan(Beta*d) )/( Z0+jZL*tan(Beta*d) ) and
% we are starting it from the end of the line and going backwards to the
% entrance, until we get to the first layer. then we will have the total
% impendace of the layers as the source "see" them. Now we can calculate
% the reflection coefficient using the formula gamma = (Zd-Z0)/(Zd+Z0) to
% find the reflection coefficient of the whole system.
% We will do the same proccess for every angle and every frequency to
% eventually draw a map of reflection coefficient for every input. 
% the lower the reflection coefficient the better the light will pass
% throuh the window and we will get clearer image.
clc
close all
clear all
N = 1:3;

c = 3e8;
%freqs for every wavelength of the visible spectrum 
f0 = (3.75e14+7.5e14)/2;
lambda0 = c/f0;
f = linspace(3.75e14,7.5e14,1000);
thetaI = linspace(-pi/4,pi/4,1000); %angle of incidance
%frequency, in our case its the center of the visible spectrum which is
%600nm
%We assume that the layers arent magntic so mu==mu -> mur=1
n = [1 zeros(1,length(N)) sqrt(2.25)]; %the first refracting index is 1 (air)
EpsilonsR = [zeros(1,length(N)) 2.25]; %EpsilonR placeholder matrix
%Epsilon & refracting index decleration
if length(N) >= 2
    n(2) = sqrt(1.257);
    n(3) = sqrt(1.773);
    EpsilonsR(1) = 1.257;
    EpsilonsR(2) = 1.773;
    if length(N) >= 3
        n(2) = sqrt(1.131);
        n(3) = sqrt(1.493);
        n(4) = sqrt(1.970);
        EpsilonsR(1) = 1.131;
        EpsilonsR(2) = 1.493;
        EpsilonsR(3) = 1.970;
        if length(N) >= 4
            n(2) = sqrt(1.0682);
            n(3) = sqrt(1.301);
            n(4) = sqrt(1.710);
            n(5) = sqrt(2.085);
            EpsilonsR(1) = 1.0682;
            EpsilonsR(2) = 1.301;
            EpsilonsR(3) = 1.710;
            EpsilonsR(4) = 2.085;
        end
    end
end
lambda = (lambda0./n);
d = lambda./4; %length of every section is quarter wavelength

for i = 1:length(f)
    for j = 1:length(thetaI)
        Zmat(1) = (120*pi)/cos(thetaI(j)); %Z of air
        dnew(1) = d(1)*cos(thetaI(j)); % d of air
        temp = thetaI(j);
        for k = 1:length(N)+1
            Zmat(k+1) = (120*pi/n(k+1))/cos(asin(n(k)*sin(temp)/n(k+1)));
            temp = asin(n(k)*sin(temp)/n(k+1));
            dnew(k+1) = d(k+1)*cos(temp);
        end
        beta(1) = 2*pi*f(i)/c;
        for k = 1:length(N)
            beta(k+1) = 2*pi*f(i)*n(k+1)/c; 
        end 
        %delete of irrelevant things
        beta(:,1) = []; 
        Za = Zmat(:,1);
        Zmat(:,1) = [];
        dnew(:,1) = [];
        for k = length(N)+1:-1:2 %looping for each elemnt in Zmat
            if k == length(N)+1 % runs on the first loop, place Z of glass
               Zd(k) = Zmat(k-1)*((Zmat(k)+1j*Zmat(k-1)*tan(dnew(k-1)*beta(k-1)))/( Zmat(k-1)+1j*Zmat(k)*tan(dnew(k-1)*beta(k-1))));
            else 
               Zd(k) = Zmat(k-1)*((Zd(k+1)+1j*Zmat(k-1)*tan(dnew(k-1)*beta(k-1)))/( Zmat(k-1)+1j*Zd(k+1)*tan(dnew(k-1)*beta(k-1))));
            end
        end   
        Zd = Zd(:,2);
        gamma(j,i) = (Zd-Za)./(Zd+Za);
    end
end
figure();
imagesc(f,flipud(thetaI),abs(gamma)); 
title('Part 4.4 N layers - TE');
ylabel(' angle [\theta]');
xlabel('Freq [Hz]');
set(gca,'YDir','normal')
colormap('jet');
colorbar('fontsize',10);