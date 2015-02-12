% Manel Soria - ETSEIAT
% Demo of aga function
% In this example, we look for a global minimum of a funcion with many local minima


function test_aga_2

% Manel Soria - ETSEIAT
% Demo of aga function
% knapsack problem solved with GA
% In this example the search space is a a set of 
% [0,1] integer values.


rng(5) % use a known random seed for repetabiliy of the problem

N = 30; % Number of items
maxval=10; % maximum possible value of one item
weights = randi(10,N,1); % Weights
values = randi(maxval,N,1); % Values
capacity = floor(sum(weights)/2);
rng('shuffle'); % we don't want repetability in the GA 

% each item can be inside the sack (1) or outside (0)
% if it is inside, contributes to value but also weight
% we look for maximum value, for total weight <= capacity

% an individual will be a vector of 0,1 values, represented using
% floating point numbers for simplicity, but aga would allow different
% data types

HUGE=N*maxval*100;

% To minimize, we change the sign
% To discard non-feasible solutions, we add a penalty
fitness=@(x) -dot(x,values)+HUGE*N*(dot(x,weights)>capacity); 
ranfun=@()  randi([0,1],[1,N]);




ng=1000; % number of generations
np=100; % population size, NOT directly used by aga but useful here
nm=2; % elite size
nr=floor(np*0.2); % number of mutants
nn=floor(np*0.9); % number of normal sons (the remaining are newcommers)
goal=-Inf; % never stop



% Now we construct an initial population list
pop={};
for i=1:np
    pop{i,1}=ranfun(); 
end
pop{np+1,1}=pop{1};



    function printpop(pop) 
        for ii=1:size(pop,1)
            pop{ii,1}
        end
    end

function popu=funique(pop) % discard identical individuals
    matu=unique(cell2mat(pop),'rows');
    popu=mat2cell(matu,ones(1,size(matu,1)));
end



function x=mutfun(x,f) 
    p=randi(size(x,2));
    x(p)=1-x(p);
end
                

x=ranfun();
y=ranfun();

function c=mate(a,b)
    e=size(a,2);
    p=randi(size(a,2));
    c=[ a(1:p-1) b(p:e) ];
end


    function prinfun(x)
        for ii=1:size(x,2)
            fprintf('%d ',x(1,ii));
        end
    end


info=10; % Verbosity level
label=10000; % Just a label in case output is to be filtered


[pop2,best,nite]=aga(info,label, ...
                 pop,...
                 ng,nm,nr,nn,goal,...
                 @funique,fitness,@mutfun,@mate,ranfun,@prinfun);
  

pop2{1}
best
nite

end
