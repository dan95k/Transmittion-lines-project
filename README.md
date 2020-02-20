Task 2 is focused around the binomial and chebychev trnasformers, using different methods of matching in order
to match an input wave to the load. In the binomial part we match using the binomial formula according to POZAR book
and fior the chebychev we used the chebychev polinomial to solve the problem

Task 3 is another approach to the same problem of matching th loads in a system, this time we used
an exponential approach which saying it is possible to match the impendaces by increasing the length
of the transmittion line.

Task 4 (Main code)
In this code we perform an analysis of a light passing through an arrayof layers in a window. 
We shoot a ray of light into the window and we want it to pass the window without it breaking. 
the ranges of the incidence angle and wavelength as follows: from -45 degrees to 45 degrees and for
the wavelength are the whole visible spectrum. 
in order to do that we use an anlogy to a transmittion line when we want to match the input
impadence to the load (input impadence is the space before the window and the load is the space after the window) 
when every layer we add is another transmittion line. Epsilons of every metter is given.

The idea is as follows, we calculate each layer zeta using the epsilons given using the forumla sqrt(mue0*muer/epsilon0*epsilonr).
Later we find each layer width as the width changes for every wavelngth and angle of incidance. 
Beta is calculated with the forumla 2*pi*f(i)*n(k)/c when f is the frequency and n is the refraction index of the layer and c
is the speed of light.
As in transmittion lines we have an characteristic impendace so we have here, 
we calculate it here (zeta of the metter)/cos(the angle in the layer)
after we have the characteristic impendaces beta and the width of every layer we can procced to calculate the reflection coefficients.
We use the formula Z0*( ZL+jZ0*tan(Beta*d) )/( Z0+jZL*tan(Beta*d) ) and we are starting it from the end of the line and going backwards to the
entrance, until we get to the first layer. then we will have the total impendace of the layers as the source "see" them. 
Now we can calculate the reflection coefficient using the formula gamma = (Zd-Z0)/(Zd+Z0) to
find the reflection coefficient of the whole system. We will do the same proccess for every angle and every frequency to
eventually draw a map of reflection coefficient for every input. 
the lower the reflection coefficient the better the light will pass throuh the window and we will get clearer image.
