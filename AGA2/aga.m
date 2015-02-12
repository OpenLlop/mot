function [lop,fop,nite,history]=aga(ninfo,label, ... 
                           pop, ... 
                           ng,N, goal, ... 
                           funique,fitfun,mutfun,reproduccio,ranfun,prifun)  

                       % 1-parfor que sigui opcional
                       % 2-si alguna funcio es empty que no la cridi si no es
                       % imprescindible (ie, print, mutacio)
                       % 3-abans de calcular, mirar si ja hem calculat
                       %
                       %-fem servir una funcio isequivalent(a,b) de l'usuari
                       % si torna 1 -> no es recalcula
                       % si torna 0 -> si es recalcula
                       % -si isequivalent es empty, internament fa servir isequuak
                       % va construint una llista de individuos coneguts, fins arribar a NCACHE
                       % si NCACHE=0, no fa res d'aixo (en algun cas sera lo millor ja que el 
                       % cost de buscarlo es creixent amn NCACHE*NP
                       
                       
% Iterates to find mimumum of a function using Genetic Algorithm v1.01
% (c) 2013 - Manel Soria - ETSEIAT 
%
% ninfo:    iteration control; prints every ninfo iterations 
% label:    integer number that precedes the prints in case output is to be
%           filtered
% pop:      list with initial population elements
% ng:       number of generations
% N:        population control parameters
%   N(1)    ne: number of elite individuals that remain unchanged
%   N(2)    nm: number of mutants
%   N(3)    nn: number of newcomers
%               the rest are descendants
%   N(4)    na: number of parents. The descendants are choosen 
%                   among the na best individuals 
%           The rest of individuals (ie: nn=length(pop)-ne+nm+nd) are
%                   newcommers, randomly choosen
% goal:     If function value is below goal, iterations are stopped
% 
% If there are less than nm-1 non-identical indivials, population is 
% considered degenerate and iterations stop
%
% Call back functions to be provided by user:
% funique:  Deletes repeated individuals in a population
%           Receives a population and returns a population
%           (a population is a list of individuals)
% fitfun:   Fitness function, given one individual returns its fitness
%           (RECALL that in this GA algoritm fitness is MINIMIZED)
% mutfun:   Mutation funcion, given one individual and its fitness,
%           mutfun should return a mutant individual. Fitness is given
%           in case mutation intensity is to be decreased when close to
%           the goal
% reproduccio: Given two individuals, returns a descendant 
% ranfun:   Returns a random individual
% prifun:   Prints individual
%
% aga returns:
% lop:      list with the population sorted by  fitness  
% fop:      minimum value of fitfun found
% nite:     number of iterations performed 
% history:  vector with the best value found after each iteration
       

if isnumeric(pop) 
    NI=pop;
    pop = cell(1,NI);
    for i=1:NI
        pop{i}=ranfun();
    end
end
    
ne=N(1); 
nm=N(2);
nn=N(3);
na=N(4);
ps=length(pop);

history=[];

nd=ps-N(1)-N(2)-N(3);
if nn<0
    error('aga: nn must be positive');
end

g=1; % generation

    while true
        nite=g;

        pop=funique(pop);
        
        if length(pop)<na 
                if (info>0)
                    fprintf('GA label=%d degenerate population\n',label);
                end
                break
        end
        
        for i=length(pop)+1:ps % repopulation
            pop{end+1}=ranfun();
        end
            
        parfor i=1:ps % parfor
            fi(i)=feval(fitfun,pop{i});
        end
        
        
        [fi,i]=sort(fi);
        pop=pop(i);

        history(end+1)=fi(1);

        if ninfo>2
            fprintf('GA label=%d g=%3d ng=%d best=%e ',label,g,ng,fi(1));
            if ~isempty(prifun)
                prifun(pop{1});
                fprintf('\n');
            else
                fprintf('\n');
            end
        
        end
        
        
        if fi(1)<=goal | g==ng 
            lop=pop;
            fop=fi(1);
            if ninfo>0 
                fprintf('GA label=%d best=%e ',label,fop);
                if ~isempty(prifun)
                    prifun(pop{1});
                    fprintf(' ');
                end             
                if fop<goal
                    fprintf('goal=%e achieved !!\n',goal);
                else
                    fprintf('max. iterations reached, leaving\n');
                end
            end
            return; %%%%%
        end
        
        % next generation:
        pos=1;
        
        for i=1:ne % elite
            pop2{pos}=pop{pos};
            pos=pos+1;
        end
        
        for i=1:nm % mutants
            pop2{pos}=mutfun(pop{pos},fi(pos));
            pos=pos+1;
        end
        
        for i=1:nd % descendants
            pare=randi([1,na]); % parents are choosen among np best 
            mare=randi([1,na]);            
            pop2{pos}=reproduccio(pop{pare},pop{mare});
            pos=pos+1;
        end       
                
        for i=1:nn % newcommers
            pop2{pos}=ranfun();
            pos=pos+1;
        end

        pop=pop2;
        
        g=g+1; % next generation
    end


    
end