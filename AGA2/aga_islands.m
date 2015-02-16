function [pops,fops,best,fbest]=aga_islands(ninfo,label,...
    pops, ngg, nemi, ng, N, goal,... %  
    funique,fitfun,mutfun,reproduccio,ranfun,prifun)

% pops: either a
%  -list of lists with population OR
%  -a two components vector
%    with the number of islands and the number of individuals per island

% lpops: last populations
% fops: fitness of the best in each island
% best: best invidual
% fbest: its fitness

if (~iscell(pops))  % Only population size is given                                   
    ni=pops(1);
    np=pops(2);
    pops={};
    for island=1:ni
        for i=1:np
            pops{island}{i}=ranfun(); 
        end
    end
else % Initial population is given .. we check it
    
    ss=size(pops);
    if (ss(1)~=1) 
        error('Wrong islands initial population');
    end
    ni=ss(2); % number of islands
    for i=1:ni
        ss=size(pops{i});
        if (ss(1)~=1) 
            error('Wrong islands initial population');
        end
        if i==1
            np=ss(2);
        else
            if np~=ss(2)
                error('Wrong population size');
            end
        end
    end
end


    
    
if ninfo>0
      fprintf('aga_islands begin ni=%d np=%d ngg=%d ng=%d \n',...
          ni,np,ngg,ng);
end

gg=1; % global generation
while true
    fops=zeros(ni,1); 
    for island=1:ni % each island evolves separately
        [lop,fops(island),nite,history]=aga(  ninfo-1,label+gg*1000+island,...
                                    pops{island},...
                                    ng,N, goal,...
                                    funique,fitfun,mutfun,reproduccio,ranfun,prifun);
        pops{island}=lop; % population at the end of local iterations
        if ninfo>1
            fprintf('aga_islands label= %d end ite global= %d illa= %d fbest= %e ',...
                label,gg,island,fops(island));
            if ~isempty(prifun) 
                fprintf('best: ');
                prifun(lop{1});
            end
            fprintf('\n');            
        end
        
    end
    
    [fbest,ibest]=min(fops); % optimal of all the islands, island where it lives
    if ninfo>0
        fprintf('aga_islands label= %d gg= %d fbest= %8.3e island= %d ',...
            label,gg,fbest,ibest);
        if ~isempty(prifun) 
            fprintf('best: ');
            prifun(lop{1});
        end
        fprintf('\n');         
    end
    
    if fbest<=goal | gg==ngg
        if ninfo>1
            fprintf('aga_islands stoping because ');
            if (gg==ngg) 
                fprintf('maximum number of global iterations %d has been reached \n', ng);
            else
                fprintf('goal %e has been reached \n', goal);
            end
        end
        best=pops{ibest}{1};
        break;
    end
    
    for island=1:ni % emigration
            dest=island;
            orig=island+1;
            if orig>ni
                orig=1;
            end
            for q=1:nemi
                mor=length(pops{island})-q+1;
                emigra=q;
                if ninfo>10
                    fprintf('aga_islands emigrates from island=%d indi=%d to replace island=%d indi=%d \n',...
                        orig,emigra,dest,mor);
                    pops{dest}{mor}=pops{orig}{emigra};       
                end
            end
    end 

    gg=gg+1;


end % global iteration