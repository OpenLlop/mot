function test_aga_arnau
% Exemple: els individuals son una estructura, i el problema d'optimizacio
% te un integer i un vector double

% Busquem una de dues coses 
% Si Q=0, un vector r3 tq norma(r)=1 i r1 = r2
% Si Q=1, un vector r3 tq norma(r)=1 i r1 = -r2

% bitxo te dos camps:
% Q
% v - vector

    function g=goal(bitxo)  % en buscarem un minim; en aquest cas pren valor 0 a l'objectiu, 
        r=bitxo.v;
        g1=norm(r)-1;
        if bitxo.Q==1 
            g2=r(1)-r(2);
        else
            g2=r(1)+r(2);
        end
        g=abs(g1)+abs(g2);
    end

    function bitxo=random_individual() % nosaltres triem que cada individuo sigui un vector fila
        min=-10;
        max=10;
        bitxo.v=min+(max-min).*rand(1,3);
        if (rand>0.5)
            bitxo.Q=1;
        else
            bitxo.Q=0;
        end
    end

    function bitxo=reproduccio(pare,mare)
        bitxo.Q=pare.Q;
        bitxo.v=(pare.v+mare.v)*0.5;
    end

    function bitxo2=mutacio(bitxo,f)
        bitxo2=bitxo;
        if (rand<0.25) 
            bitxo2.Q=1-bitxo.Q;
        else
            bitxo2.v=bitxo.v+(rand(1,3)-rand(1,3))*1e-3;
        end
    end

    funique=@(pop) pop; % Given a population, returns a population
                        % unique should discard identical individuals, 
                        % Here, we ignore this possibility as the individuals are R^2 vectors

    prifun=@(x) fprintf('%d % e % e % e',x.Q,x.v(1),x.v(2),x.v(3));


np=1000;  % population size (not directly used by aga)

N=[1,... % elite: doesn't change
   floor(np*0.1),... % number of mutants
   floor(np*0.05),...% number of newcommers
   floor(np*0.2) % number of parents 
    ];

ng=50; % number of generations


stop=1e-5; % goal (if aga finds a value below goal, iterations will be stopped)


rng('shuffle'); % IMPORTANT ! we don't want repetability in the GA 


info=5; % Verbosity level
label=10000; % Just a label (an arbitray integer number) 
             % in case output is to be filtered

[pop2,best,nite,history]=aga(info,label, ...
                 np,... % initial population size; if a list is given, it is used as initial population
                 ng,N,stop,...
                 funique,@goal,@mutacio,@reproduccio,@random_individual,prifun);
    % pop2 is the sorted population after iterating; pop2{1} is the best
    %   solution available
    % best=fitfun(pop2{1}) is the optimal value
    % nite is the number of iterations performed
found=pop2{1};

figure, semilogy(history,'o-');

found
best

% Now, we can easily improve the accuracy of the local extremum found
%found_fminsearch=fminsearch(fitfun,found)


end


