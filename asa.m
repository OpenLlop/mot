function [ bestind, bestfit, nite, history ] = asa ( ninfo, label, ...
    A0, nitemax, mu, goal, ...
    fitfun, mutfun, prifun )
% Iterates to find mimumum of a function using Simulated Annealing
% (c) 2013 - Manel Soria - ETSEIAT - v1.0
% (c) 2015 - Manel Soria, David de la Torre - ETSEIAT - v1.01
%
% ninfo:    prints information every ninfo iterations
% label:    integer number preceding the prints in case output is to be
%           filtered
% A0:       initial guess
% nitemax:  maxim number of iteracions
% mu:       Simulated annealing parameter (eg, 0.2, read below)
% goal:     If function value is below goal, iterations are stopped
%
% Call back functions to be provided by user:
% fitfun:   fitness function
% mutfun:   mutation (change) function. Receives an individual and its
%           fitness and returns a modified individual
% prifun:   Prints individual
%
% Returns:
% bestind:  Best individual
% bestfir:  Fitness of the best individual
%
% About mu parameter (from
% http://sourceforge.net/p/sbsi/discussion/1048774/thread/6caecc11/)
% High values of Mu mean that the probability of accepting a worse solution
% from one iteration to the next is low, and hence the solution is likely 
% to head straight to a local minimum. Conversely, low values of Mu mean 
% that the probability of accepting a worse solution from one iteration to 
% the next is high, and hence the solution is able to jump out of local 
% minima but may also not be able to converge at all. So finding a good 
% value of Mu is a tradeoff between these extremes and the absolute value 
% needs to be empirically determined for different models.

% Print extra info (for debugging purposes)
extrainfo = 0;

% History
history = [];

% Preprocessing
A = A0; % Set initial guess
fita = fitfun(A); % Compute initial guess fitness
bestind = A; % Best individual
bestfit = fita; % Best fitness

% Iterate until convergence or max iterations
for nite=1:nitemax
    
    % Must print something
    mustprint = (nite==1 || mod(nite,ninfo)==0);
    
    % Print info
    if mustprint
        fprintf('SA label=%d nite=%d ',label,nite); % Print info
        prifun(A); % Print individual
        fprintf(' fitA=%f fitbest=%f \n',fita,bestfit); % Print fitness
    end;

    % Break iteration if convergence is achieved
    if fita<goal
        break;
    end;
    
    % Mutate new individual B
    B = mutfun(A,fita);
    
    % Compute B fitness
    fitb = fitfun(B);

    % Save best individual/fitness
    if fitb<bestfit
        bestind = B; % Best individual
        bestfit = fitb; % Best fitness
    end;
    
    % Save fitness history
    history(nite) = bestfit; %#ok

    % Compute fitness difference between B and A
    deltafit = fitb - fita;
    
    % Simulate thermal transition probability (~ exp(-Delta/(k_b*T))
    probability = exp(-deltafit / (mu * abs(fita)));
    
    % Print extra info
    if mustprint && extrainfo
        fprintf('   fitb=%8.2e deltafit=%+8.2e deltafit/fit=%+8.2e probability=%+8.2e  ',...
                    fitb,      deltafit,       deltafit/abs(fita), probability);  
    end;
    
    % Jump from A to B (randomly with thermal transition probability)
    if rand<probability % Jump
        A = B; % Jump from A to B
        fita = fitb; % Update fitness value
        if mustprint && extrainfo % Print extra info
            fprintf('jump '); 
            if deltafit<=0, fprintf('>= \n'); % Fitness improved
            else fprintf('< \n'); % Fitness worsened
            end;
        end;
    else % Do not jump
        if mustprint && extrainfo % Print extra info
            fprintf('= \n'); % Fitness does not change
        end;
    end;
    
end;

end

