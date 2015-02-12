% Manel Soria - ETSEIAT
% Demo of aga function
% In this example our seach space are regular polygons with radius <= Rmax
% and number of sizes <= Nmax
% Our goal is to find the polygon with largest surface/perimeter ratio,
% that we know it has to be: N=Nmax, R=Rmax (the most similar to the
% maximum circle

clear

Nmin=3
Nmax=18
Rmax=100

Parea=@(N,R) 0.5*R^2*N*sin(2*pi/N);
Pperimeter=@(N,R) N*2*R*sin(pi/N);

ite=@(x,a,b) (x==1)*a+(x~=1)*b; % if x then a else b


opt=@(N,R) ite(N>=Nmin && N<=Nmax && R<=Rmax,Parea(N,R)/Pperimeter(N,R),-Parea(N,R)/Pperimeter(N,R));



% Here begin ga parameters and call back functions 

ng=50; % number of generations
np=300; % population size, NOT directly used by aga but useful here
nm=2; % elite size
nr=floor(np*0.2); % number of mutants
nn=floor(np*0.9); % number of normal sons (the remaining are newcommers)
goal=-Inf; % goal (if aga finds a value below goal, iterations will be stopped)



% x(1): should be integer, number of sides
% x(2): radius
fitfun=@(x) -opt(x(1),x(2));

funique=@(pop) pop;


ranfun=@() [randi([Nmin,Nmax]),Rmax*rand];

deltai=[0,0,0,0,-1,1]; % fem que sigui mes probable deixar N igual 

mutfun=@(x,f) [x(1)+deltai(randi([1,length(deltai)])),x(2)+0.1*(rand-rand)];

reproduccio=@(x,y) ite(rand>0.5, ...
                        [x(1),0.5*(x(2)+y(2))] , ... 
                        [y(1),0.5*(x(2)+y(2))] );
                    

prifun=@(x) fprintf('%d %f',x(1),x(2));

rng('shuffle'); % we don't want repetability in the GA 


% We construct an initial population list
pop={};
for i=1:np
    pop{i}=ranfun(); 
end


info=1; % Verbosity level
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
found=pop2{1}
