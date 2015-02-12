function [popt,fopt,best,fbest]=aga_islands(ninfo,label,...
                                          niteg,ni,np, ng,nm,nr,nn, goal,...
                                          funique,fitfun,mutfun,reproduccio,ranfun,prifun)
% Iterates to find mimumum of a function using Genetic Algorithm with
% islands
% (c) 2013 - Manel Soria - ETSEIAT

% ninfo:    iteration control; prints every ninfo iterations 
% label:    integer number that precedes the prints in case output is to be
%           filtered
% niteg:    number of global iterations
% ni:       number of islands
% np:       nombre of individuals per island
%
% on each island:
% ng:       number of generations
% nm:       control of elite individuals: after sorting by fitness, 
%           pop{1} .. pop{nm-1} remain unchanged
% nr:       control of mutations: pop{nm} to pop{nr-1} become mutations 
%           of the elite 
% nn:       control of reproduction: pop{nr}..pop{nn-1} are descendants of 
%           random individuals selected from all the population
%           the remaining pop{nn}..pop{end} are newcommers, random 
%           individuals 
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
% aga_ilands returns:
% popt:     Best individual of each island
% fopt:     Fitness of the best individual of each island
% best:     Best individual overall
% fbest:    Fitness of the best individual overall

pop={};
for illa=1:ni
    for i=1:np
        pop{illa}{i}=ranfun(); 
    end
end

for g=1:niteg
    fop=zeros(ni); 
    for illa=1:ni % each island evolves separately
        [lop,fop(illa),nite]=aga(ninfo,label+g*1000+illa,pop{illa},ng,nm,nr,nn, goal,funique,fitfun,mutfun,reproduccio,ranfun,prifun);
        pop{illa}=lop;
        if ninfo<1000
            fprintf('label= %d ite global= %d illa= %d millor valor= %e ',label,g,illa,fop(illa));
            if ~isempty(prifun) 
                fprintf('best ');
                prifun(lop{1});
                fprintf('\n');
            else
                fprintf('\n');
            end
        end
        
    end
    best=min(fop);
    if best<goal
        break;
    end
    if ninfo<1000, fprintf('not solved yet .. migrations \n');
    
    for illa=1:ni % best in each island replaces the worse in the following
        if illa>1
            origen=illa-1;
        else
            origen=ni;
        end
        pop{illa}{end}=pop{origen}{1};
    end
end

popt={};
fopt=[];
for illa=1:ni
    popt{illa}=pop{illa}{1};
    fopt(illa)=fitfun(popt{illa});
end

[fopt,indx]=sort(fopt);
popt=popt(indx);

best=popt{1};
fbest=fitfun(best);

end