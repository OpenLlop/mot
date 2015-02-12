function [lop,fop,nite]=aga(ninfo,label, ... 
                           pop, ... 
                           ng,nm,nr,nn, goal, ... 
                           funique,fitfun,mutfun,reproduccio,ranfun,prifun)  
% Iterates to find mimumum of a function using Genetic Algorithm
% (c) 2013 - Manel Soria - ETSEIAT
%
% ninfo:    iteration control; prints every ninfo iterations 
% label:    integer number that precedes the prints in case output is to be
%           filtered
% pop:      list with initial population elements
% ng:       number of generations
% nm:       control of elite individuals: after sorting by fitness, 
%           pop{1} .. pop{nm-1} remain unchanged
% nr:       control of mutations: pop{nm} to pop{nr-1} become mutations 
%           of the elite 
% nn:       control of reproduction: pop{nr}..pop{nn-1} are descendants of 
%           random individuals selected from all the population
%           the remaining pop{nn}..pop{end} are newcommers, random 
%           individuals 
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
                
np=length(pop);

    for g=1:ng
        nite=g;

        pop=funique(pop);
        
        if length(pop)<nm-1
                fprintf('GA label=%d degenerate population\n',label);
                break
        end
        
        for i=length(pop)+1:np % repopulation
            pop{end+1}=ranfun();
        end
            
        parfor i=1:np
            fi(i)=feval(fitfun,pop{i});
        end
                
        [fi,i]=sort(fi);
        pop=pop(i);

        if g==1 || g==ng || mod(g,ninfo)==0
            fprintf('GA label=%d g=%3d ng=%d nm=%d nr=%d best=%e ',label,g,ng,nm,nr,fi(1));
            if ~isempty(prifun)
                prifun(pop{1});
                fprintf('\n');
            else
                fprintf('\n');
            end
        
        end
        
        
        if fi(1)<goal || g==ng 
            lop=pop;
            fop=fi(1);            
            if fi(1)<goal
                fprintf('GA label=%d goal=%e achieved !!\n',label,goal);
            else
                fprintf('GA label=%d goal=%e lograt=%e max. number of iteracions, leaving\n',label,goal,fop);
            end

            return;
        end
        

        % mutation
        for i=nm:nr-1
            mutant=randi([1,nm]);
            pop{i}=mutfun(pop{i},fi(i));
        end
        
        % reproduction
        for i=nr:nn-1 
            pare=randi([1,np]);
            mare=randi([1,np]);
            pop{i}=reproduccio(pop{pare},pop{mare});
        end
        
        % newcommers
        for i=nn:np
            pop{i}=ranfun();
        end

    end

    parfor i=1:np
        fi(i)=feval(fitfun,pop{i});
    end
    
    [fi,i]=sort(fi);
    pop=pop(i); 
 
    lop=pop; 
    fop=fi(1);
    
end