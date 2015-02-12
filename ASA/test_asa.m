% Manel Soria - ETSEIAT
% Demo of asa function
% In this example, we look for a global minimum of a funcion with many local minima

clear
close all

% Our test is a R^2->R function based on Rastrigin function. It is 
% challenging because it has infinite local extrema, located at integer
% numbers (ie, 8,-9) 
% The global minimum is at (1,1), and its value is 0

ras=@(x,y) 20+(x-1).^2+(y-1).^2-10*(cos(2*pi*(x-1))+cos(2*pi*(y-1)));


ffit=@(x) ras(x(1),x(2));

ranrange=@(a,b,n) a+(b-a)*rand(n,1); % n values uniform between a i b
neighbour=@(x,f) x+ranrange(-0.1,0.1,2); % neighbour, just a small random movement (its fitness is ignored)

initial_guess=[0.23;1.45];

prifun=@(x) fprintf('%f %f ',x(1),x(2));

[best,fit_best]=asa(100,11,initial_guess,1000,0.2,1e-3,ffit,neighbour,prifun);


best
fit_best
