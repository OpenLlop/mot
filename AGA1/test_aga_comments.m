% Manel Soria - ETSEIAT
% Demo of aga function
% In this example, we look for a global minimum of a funcion with many local minima
clear; close all;

% Our test is a R^2->R function based on Rastrigin function. It is 
% challenging because it has infinite local extrema, located at integer
% numbers (ie, 8,-9) 
% The global minimum is at (1,1), and its value is 0
ras=@(x,y) 20+(x-1).^2+(y-1).^2-10*(cos(2*pi*(x-1))+cos(2*pi*(y-1)));

% Plot the rastrigin function
% [X,Y] = meshgrid(-5:.05:5, -5:.05:5);
% surfc(-5:.05:5, -5:.05:5,ras(X,Y));
% shading interp;

% Define the GA parameters
ng = 10; % Number of generations
np = 300; % Population size, NOT directly used by aga but useful here
nm = 2; % Elite size
nr = floor(np*0.2); % Number of mutants
nn = floor(np*0.9); % Number of normal sons (the remaining are newcommers)
goal = 1e-5; % Target fitness value (if aga finds a value below goal, iterations will be stopped)

% Auxiliary function
ranrange=@(a,b,n) a+(b-a)*rand(n,1); % n values uniform between a i b

% Define GA functions
funique = @(pop) pop; % unique should discard identical individuals, we ignore this possibility as the individuals are R^2 vectors
fitfun = @(x) ras(x(1),x(2));  % fitness function - TO BE MINIMIZED. Note that its single argument is a 2-component vector
mutfun = @(x,f) x+ranrange(-0.1,0.1,2); % mutation, just a small random movement
reprofun = @(x,y) (x+y)/2; % reproduction, in this case just the average
ranfun = @() ranrange(-100,100,2); % random individual (uniform between -100 and 100)
prifun = @(x) fprintf('%f %f',x(1),x(2)); % Print an individual

% Randomize random seed
rng('shuffle'); % We don't want repetability in the GA 

% Construct an initial population list
pop = cell(np,1);
for i=1:np
    pop{i} = ranfun(); 
end;

% Define GA info data
info = 0; % Verbosity level
label = 10000; % Just a label (an arbitray integer number) in case output is to be filtered

% Execute Genetic Algorithm
[pop2,bestValGA,nite]=aga_comments(info,label, ...
                 pop,...
                 ng,nm,nr,nn,goal,...
                 funique,fitfun,mutfun,reprofun,ranfun,prifun);
    % pop2 is the sorted population after iterating; pop2{1} is the best solution available
    % best = fitfun(pop2{1}) is the optimal value
    % nite is the number of iterations performed
bestIndGA = pop2{1};

% Now, we can easily improve the accuracy of the local extremum found
options = optimset('TolFun',1e-8,'Display','none');
[bestIndFMS,bestValFMS] = fminsearch(fitfun,bestIndGA,options);

% Display values
fprintf('Algorithm \tBest individual (x,y) \tValue\n');
fprintf('GA \t\t\t%1.6f,%1.6f \t\t%1.6f\n',bestIndGA,bestValGA);
fprintf('FMS \t\t%1.6f,%1.6f \t\t%1.6f\n',bestIndFMS,bestValFMS);


