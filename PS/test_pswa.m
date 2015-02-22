% Manel Soria - ETSEIAT
% Demo of aga function
% In this example, we look for a global minimum of a funcion with 
% many local minima

clear
close all

% Our test is a R^2->R function based on Rastrigin function. It is 
% challenging because it has infinite local extrema, located at integer
% numbers (ie, 8,-9) 
% The global minimum is at (1,1), and its value is 0

ras=@(x,y) 20+(x-1).^2+(y-1).^2-10*(cos(2*pi*(x-1))+cos(2*pi*(y-1)))

% We plot it..
% [X,Y] = meshgrid(-5:.05:5, -5:.05:5);
% ras=@(x,y) 20+(x-1).^2+(y-1).^2-10*(cos(2*pi*(x-1))+cos(2*pi*(y-1)))
% surfc(-5:.05:5, -5:.05:5,ras(X,Y))
% shading interp;
% view(0,90);


np=50;  % population size (not directly used by aga)



goal=1e-5; % goal (if aga finds a value below goal, iterations will be stopped)


% auxiliary function
ranrange=@(a,b,n) a+(b-a)*rand(n,1); % n values uniform between a i b



funique=@(pop) pop; % unique should discard identical individuals, 
                    % we ignore this possibility as the individuals are R^2 vectors
fitfun=@(x) ras(x(1),x(2));  % fitness function - TO BE MINIMIZED 
                             % note that its single argument is a 2-component
                             % vector

%fitfun=@(x) x(1)^2+x(2)^2; another (easier) function

ranfun=@() ranrange(-50,50,2); % random individual (uniform between -100 and 100)
prifun=@(x) fprintf('% f  % f',x(1),x(2));

rng('shuffle'); % we don't want repetability in the GA 


% We construct an initial population list
pop={};
v={};
for i=1:np
    pop{i}=ranfun(); 
    v{i}=0.0*ranfun(); % in this algorithm, an initial (random) velocity has to 
end


% 
[best,fbest]=pswa(pop,v,9,20,5e-3,100,goal,fitfun,prifun);




