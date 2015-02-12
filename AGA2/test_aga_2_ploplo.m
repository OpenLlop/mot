function test_aga_2_ploplo

    close all
    hold off
    
    %{
    for test=1:10
        NN=6
        for qq=1:32
            val=aga_performance(NN);
            
            if val>0
                psize(qq)=NN;
                effort(qq)=val;
                bf(qq)=2^NN;
            else
                break;
            end
            
            NN=NN+1;           
        end


        semilogy(psize,effort,'o');
        if test==1
            hold on;
        end

    end
    
    semilogy(psize,bf,'r');
%}
    
    close all
    
    
    aga_performance(100)
    
    function nite_np=aga_performance(NN)
    % Manel Soria - ETSEIAT
    % Plots iterations*population size vs. problem size
    % needed by GA to solve 
    % a random knapsack problem of size NN


    rng(5) % use a known random seed for repetabiliy of the problem
    % NN = 30; % Number of items in the knapsack
    maxval=100; % maximum possible value of one item
    maxw=20; % max weight 
    weights = randi(maxw,NN,1); % Weights
    values = randi(maxval,NN,1); % Values
    capacity = floor(sum(weights)/2); % Capacity


    % we compute the solution of this problem

    [best amount] = knapsack(weights, values, capacity); 

        best
        
    rng('shuffle'); % we don't want repetability in the GA 


    % each item can be inside the sack (1) or outside (0)
    % if it is inside, contributes to value but also adds weight
    % we look for maximum value, for total weight <= capacity

    % an individual will be a vector of 0,1 values, represented using
    % floating point numbers for simplicity, but aga would allow different
    % data types
    maxval=1000;

    HUGE=NN*maxval*100;

    ng=500; % max number of generations
    np=300; % population size, NOT directly used by aga but useful here

    N=[5,... % elite
       floor(np*0.4),... % mutants
       floor(np*0.05),...% newcommers
       floor(np*0.1)
        ];
    goal=-best; 



    % Now we construct an initial population list
    pop={};
    for i=1:np
        pop{i}=ranfun(); 
    end


    info=10; % Verbosity level
    label=10000+NN; % Just a label in case output is to be filtered



    [pop2,bestga,nite,history]=aga(info,label, ...
                     pop,...
                     ng,N,goal,...
                     @funique,@fitness,@mutfun,@mate,@ranfun,@prinfun);


    plot(-history,'o');
    hold on
    plot(best*ones(length(history),1),'r');
    grid
    
    if (norm(pop2{1}-amount)==0) 
        fprintf('ok, found correct solution in nite=%d iterations\n',nite);
    else 
        nite_np=-1;

    end

    nite_np=nite*np;

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

end


