% Manel Soria - ETSEIAT
% Demo of aga_islands function
% Here we reconsider our first example but using islands

clear
close all

% These three lines are the main difference from the test_aga code:
ni=20; % number of islands
ngg=10; % number of global iterations
np=40; % population of each island de cada illa
ng=10; % number of local generations
nemi=5; % number of emigrants

% All this is more or less like in the single island version:

label=10000; % arbitrary integer to filter solution

% Our test is a R^2->R function based on Rastrigin function. It is 
% challenging because it has infinite local extrema, located at integer
% numbers (ie, 8,-9) 
% The global minimum is at (1,1), and its value is 0

ras=@(x,y) 20+(x-1).^2+(y-1).^2-10*(cos(2*pi*(x-1))+cos(2*pi*(y-1)))

% We plot it..
%{
[X,Y] = meshgrid(-5:.05:5, -5:.05:5);
ras=@(x,y) 20+(x-1).^2+(y-1).^2-10*(cos(2*pi*(x-1))+cos(2*pi*(y-1)))
surfc(-5:.05:5, -5:.05:5,ras(X,Y))
shading interp;
%}

N=[1,... % elite
   floor(np*0.1),... % mutants
   floor(np*0.1),...% newcommers
   floor(np*0.2)
    ];

niteg=3; % number of global iterations
ng=5; % number of generations
goal=1e-5; % goal (if it finds a value below goal, iterations will be stopped)

% auxiliary function
ranrange=@(a,b,n) a+(b-a)*rand(n,1); % n values uniform between a i b


funique=@(pop) pop;
fitfun=@(x) ras(x(1),x(2)); 
mutfun=@(x,f) x+ranrange(-0.1,0.1,2);
reproduccio=@(x,y) (x+y)/2;
ranfun=@() ranrange(-100,100,2);
prifun=@(x) fprintf('%8.3f %8.3f',x(1),x(2));


% We can just give the number of islands and individuals,
% aga_islands generates the populations calling our ranfun
% pops=[ni np]; 

% or we can generate our own initial populations, for instance to begin with
% a previous computation:

pops={};
for illa=1:ni
    for i=1:np
        pops{illa}{i}=ranfun(); 
    end
end
    
[popt,fopt,best,fbest]=aga_islands(1,label,... % verbosity parameters
                pops, ...
                ngg, nemi,ng,N,goal,...
                funique,fitfun,mutfun,reproduccio,ranfun,prifun); % callbacks

best
fbest

%{
for illa=1:ni
    for ind=1:np
        fprintf('island=%2d ind=%2d ',illa,ind); 
        prifun(popt{illa}{ind});
        fprintf(' fit=%e \n',fitfun(popt{illa}{ind}));
    end
end
%}
% partint d'un punt proper a l'optim, amb un metode convencional la clavem
%options=optimset('TolFun',1e-6,'Display','iter');
%sol=fminsearch(fitfun,best,options)


