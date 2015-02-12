% Manel Soria - ETSEIAT
% Demo of aga function
% In this example, we look for a global minimum of a funcion with many local minima

clear
close all

% Our test is a R^2->R function based on Rastrigin function. It is 
% challenging because it has infinite local extrema, located at integer
% numbers (ie, 8,-9) 
% The global minimum is at (1,1), and its value is 0

ras=@(x,y) 20+(x-1).^2+(y-1).^2-10*(cos(2*pi*(x-1))+cos(2*pi*(y-1)))

% We plot it..
%
[X,Y] = meshgrid(-5:.05:5, -5:.05:5);
ras=@(x,y) 20+(x-1).^2+(y-1).^2-10*(cos(2*pi*(x-1))+cos(2*pi*(y-1)))
surfc(-5:.05:5, -5:.05:5,ras(X,Y))
shading interp;
%

np=200;  % population size (not directly used by aga)

N=[1,... % elite: doesn't change
   floor(np*0.1),... % number of mutants
   floor(np*0.05),...% number of newcommers
   floor(np*0.2) % number of parents 
    ];

ng=50; % number of generations


goal=1e-5; % goal (if aga finds a value below goal, iterations will be stopped)


% auxiliary function
ranrange=@(a,b,n) a+(b-a)*rand(n,1); % n values uniform between a i b



funique=@(pop) pop; % Given a population, returns a population
                    % unique should discard identical individuals, 
                    % Here, we ignore this possibility as the individuals are R^2 vectors
                    
fitfun=@(x) ras(x(1),x(2));  % fitness function - TO BE MINIMIZED 
                             % note that its single argument is a 2-component
                             % vector
                             
mutfun=@(x,f) x+ranrange(-0.1,0.1,2); % mutation, just a small random movement
reproduccio=@(x,y) 0.9*x+0.1*y; % reproduction, just average
ranfun=@() ranrange(-200,200,2); % random individual (uniform between -100 and 100)
prifun=@(x) fprintf('% f  % f',x(1),x(2));

rng('shuffle'); % we don't want repetability in the GA 


info=5; % Verbosity level
label=10000; % Just a label (an arbitray integer number) 
             % in case output is to be filtered

[pop2,best,nite,history]=aga(info,label, ...
                 np,... % initial population size; if a list is given, it is used as initial population
                 ng,N,goal,...
                 funique,fitfun,mutfun,reproduccio,ranfun,prifun);
    % pop2 is the sorted population after iterating; pop2{1} is the best
    %   solution available
    % best=fitfun(pop2{1}) is the optimal value
    % nite is the number of iterations performed
found=pop2{1};

figure, semilogy(history,'o-');

found
best

% Now, we can easily improve the accuracy of the local extremum found
found_fminsearch=fminsearch(fitfun,found)



