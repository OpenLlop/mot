function test_aga_2
% Manel Soria - ETSEIAT
% Demo of aga function
% knapsack problem solved with GA
% In this example the search space is a a set of 
% [0,1] integer values.


rng(5) % use a known random seed for repetabiliy of the problem
NN = 64; % Number of items in the knapsack
maxval=100; % maximum possible value of one item
maxw=20; % max weight 
weights = randi(maxw,NN,1); % Weights
values = randi(maxval,NN,1); % Values
capacity = floor(sum(weights)/2); % Capacity

close all

rng('shuffle'); % we don't want repetability in the GA 


% each item can be inside the sack (1) or outside (0)
% if it is inside, contributes to value but also adds weight
% we look for maximum value, for total weight <= capacity

% an individual will be a vector of 0,1 values, represented using
% floating point numbers for simplicity, but aga would allow different
% data types
maxval=1000;

HUGE=NN*maxval*100;

ng=300; % number of generations
np=100; % population size, NOT directly used by aga but useful here

N=[5,... % elite
   floor(np*0.4),... % mutants
   floor(np*0.05),...% newcommers
   floor(np*0.1)
    ];
goal=-Inf; 



% Now we construct an initial population list
pop={};
for i=1:np
    pop{i}=ranfun(); 
end


info=10; % Verbosity level
label=10000; % Just a label in case output is to be filtered



[pop2,best,nite,history]=aga(info,label, ...
                 pop,...
                 ng,N,goal,...
                 @funique,@fitness,@mutfun,@mate,@ranfun,@prinfun);
  

prinfun(pop2{1});



plot(-history,'o-')


    function r=ranfun() 
        r=transpose(randi([0,1],[1,NN]));
    end

    function f=fitness(x) 
        f=-dot(x,values)+HUGE*NN*(dot(x,weights)>capacity);
    end

    function printpop(pop) 
        for ii=1:length(pop)
            for j=1:length(pop{ii})
                fprintf('%d',pop{ii}(j));
            end
            fprintf('\n');
        end
    end

    function popu=funique(pop) % discard identical individuals
        mat=transpose(reshape(cell2mat(pop),NN,np)); % convert population to a matrix, individuals by rows
        matu=transpose(unique(mat,'rows'));
        popu=mat2cell(matu,NN,ones(1,size(matu,2)));
    end

    function x=mutfun(x,f) 
        for i=1:4
            p=randi(size(x,1));
            x(p)=1-x(p);
        end
    end
                

    function c=mate(a,b)
        e=size(a,1);
        p=randi(size(a,1));
        c=[ a(1:p-1); b(p:e) ];
    end


    function prinfun(x)
        for ii=1:size(x,1)
            fprintf('%d',x(ii,1));
        end
        fprintf(' ');
    end
end
