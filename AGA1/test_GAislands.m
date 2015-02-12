% Manel Soria - ETSEIAT
% Demo of aga_islands function
% Here we reconsider our first example but using islands

clear
close all

% These three lines are the main difference from the test_aga code:
nilles=4; % number of islands
niteglob=3; % number of global iterations
np=100; % population of each island de cada illa

% All this is more or less like in the single island version:

label=10000; % arbitrary integer to filter solution

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

% parameters for each island
ng=50; % number of generations in each island
nm=2; % elite size
nr=floor(np*0.2); % number of mutants
nn=floor(np*0.9); % number of normal sons (the remaining are newcommers)

goal=1e-5; % goal (if it finds a value below goal, iterations will be stopped)

% auxiliary function
ranrange=@(a,b,n) a+(b-a)*rand(n,1); % n values uniform between a i b


info=1000

funique=@(pop) pop;
fitfun=@(x) ras(x(1),x(2)); 
mutfun=@(x,f) x+ranrange(-0.1,0.1,2);
reproduccio=@(x,y) (x+y)/2;
ranfun=@() ranrange(-100,100,2);
prifun=@(x) fprintf('%f %f',x(1),x(2));

% Here we don't need to create an initial population
[popt,fopt,best,fbest]=aga_islands(1000,label,... % verbosity parameters
                niteglob,nilles,np,  ng,nm,nr,nn,  goal,  ... 
                funique,fitfun,mutfun,reproduccio,ranfun,prifun); % callbacks

best

% partint d'un punt proper a l'optim, amb un metode convencional la clavem
%options=optimset('TolFun',1e-6,'Display','iter');
%sol=fminsearch(fitfun,best,options)


