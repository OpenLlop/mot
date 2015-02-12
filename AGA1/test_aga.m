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
[X,Y] = meshgrid(-5:.05:5, -5:.05:5);
ras=@(x,y) 20+(x-1).^2+(y-1).^2-10*(cos(2*pi*(x-1))+cos(2*pi*(y-1)))
surfc(-5:.05:5, -5:.05:5,ras(X,Y))
shading interp;


ng=10; % number of generations
np=300; % population size, NOT directly used by aga but useful here
nm=2; % elite size
nr=floor(np*0.2); % number of mutants
nn=floor(np*0.9); % number of normal sons (the remaining are newcommers)
goal=1e-5; % goal (if aga finds a value below goal, iterations will be stopped)


% auxiliary function
ranrange=@(a,b,n) a+(b-a)*rand(n,1); % n values uniform between a i b



funique=@(pop) pop; % unique should discard identical individuals, 
                    % we ignore this possibility as the individuals are R^2 vectors
fitfun=@(x) ras(x(1),x(2));  % fitness function - TO BE MINIMIZED 
                             % note that its single argument is a 2-component
                             % vector
mutfun=@(x,f) x+ranrange(-0.1,0.1,2); % mutation, just a small random movement
reproduccio=@(x,y) (x+y)/2; % reproduction, just average
ranfun=@() ranrange(-100,100,2); % random individual (uniform between -100 and 100)
prifun=@(x) fprintf('%f %f',x(1),x(2));

rng('shuffle'); % we don't want repetability in the GA 


% We construct an initial population list
pop={};
for i=1:np
    pop{i}=ranfun(); 
end


info=10; % Verbosity level
label=10000; % Just a label (an arbitray integer number) 
             % in case output is to be filtered

[pop2,best,nite]=aga(info,label, ...
                 pop,...
                 ng,nm,nr,nn,goal,...
                 funique,fitfun,mutfun,reproduccio,ranfun,prifun);
    % pop2 is the sorted population after iterating; pop2{1} is the best
    %   solution available
    % best=fitfun(pop2{1}) is the optimal value
    % nite is the number of iterations performed
found=pop2{1};


% Solution with simulated annealing
%
%             [found,best]=asa(10000,label,...
%                 ranfun(),1000000,0.2,goal,fitfun,mutfun,prifun);


found
best

% Now, we can easily improve the accuracy of the local extremum found
options=optimset('TolFun',1e-8,'Display','none');
found_fminsearch=fminsearch(fitfun,found,options)



