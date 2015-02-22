function [ best, fbest ] = pswa ( pop, v, c1, c2, P1, ... 
                            nite, goal, ... 
                            fitfun, prifun)
% Iterates to find minimum of a function using Genetic Algorithm
% (c) 2015 - Manel Soria, David de la Torre - ETSEIAT - v1.00
%
% c1: local learning factor (about 2)
% c2: global learning factor (about 2)
% P1: maximum step size (velocity) allowed in one iteration
%     expressed as a fraction of the domain size
%     The domain size is estimated based on the initial
%     distribution of the particles
%     P1 should be about 1e-2
%     Low values of P1 slow down the algorithm, too high values can be
%     unstable
% a ranfun should be given to generate initial population if empty
% nite, goal, as in GA
% es podria tornar lo mateix que en GA, o sigui la poblacio final, evoluci?
% del residuo etc per fer-les lo m?s similars possible

ps=length(pop);
nc=length(pop{1}); % population elements are (must be) assumed to be vectors
fi=zeros(ps,1);   
fib=zeros(ps,1); % personal best of each particle
pbest=pop;
ite=1;

best=pop{1}; % we begin with an arbitrary best
fbest=0; 

% Estimate the domain size to limit particle velocity
MAX=pop{1}; % arbitrary
MIN=pop{1};
for i=2:ps
    MAX=max(MAX,pop{i});
    MIN=min(MIN,pop{i});
end
MAXV=P1*(norm(MAX-MIN)); % maximum velocity 

while true

    for i=1:ps % evaluate fitness
        fi(i)=feval(fitfun,pop{i});
        if ite==1 || fi(i)<fib(i) % update personal best of each particle
            fib(i)=fi(i);
            pbest{i}=pop{i};
        end
    end
    
    
    [m,pm]=min(fi); 

    if (m<fbest || ite==1) 
        best=pop{pm};
        fbest=m;
    end
    
    % this should be improved, "info" implemented etc
    % however note that prifun can also print velocities as they
    % are of the same type as individuals
    fprintf('ite=%d best=',ite); 
    prifun(pop{pm});
    fprintf(' vbest= ');
    prifun(v{pm});
    fprintf(' fbest=%e \n',m);
    
    
    if m<goal 
        best=pop{pm};
        fbest=m;
        break
    end
    
    if ite>nite
        % ok=0 % set flag (not implemented !!)
        break;
    end
    
    for i=1:ps % update positions
        pop{i}=pop{i}+v{i};
    end

    for i=1:ps % update velocities
        v{i}=v{i}+... % previous velocity
             c1*rand*(pbest{i}-pop{i})+... % local learning 
             c2*rand*(best-pop{i}); % global learning
    end
   

    for i=1:ps % limit velocities
        nvi=norm(v{i});
        if nvi>MAXV
            v{i}=MAXV*v{i}/nvi;
        end
    end
    
    % 
    %plot .. commented out .. perhaps a plotpopulation function could be
    %provided by the user ?
    hold off;
    for i=1:ps
        q=pop{i};
        plot(q(1),q(2),'o');
        hold on;       
    end
    axis([-10 10 -10 10]);
    drawnow;
    %
   
    ite=ite+1;
end


end