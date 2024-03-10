function [ bestind, bestfit, nite, lastpop, lastfit, history ] = asa ( ...
    opts, A0, goal, nitemax, mu, fitfun, mutfun, prifun )
%ASA finds mimumum of a function using Simulated Annealing (SA)
%
%Programmers:   Manel Soria         (UPC/ETSEIAT)
%               David de la Torre   (UPC/ETSEIAT)
%               Arnau Miro          (UPC/ETSEIAT)
%Date:          10/05/2018
%Revision:      3
%
%Usage:         [bestind, bestfit, nite, lastpop, lastfit, history] = ...
%                   ASA ( opts, A0, goal, DATA )
%
%Inputs:
%   opts:       function control parameters [struct] (optional)
%       ninfo:  prints information (0=none, 1=basic, 2=extended)
%       label:  integer number that precedes the prints in case output is
%               to be filtered
%       fhist:  saved history level (0=none, 1=just fitness, 2=all data)
%       	  0:      history = []
%         	1:      history(nite) = bestfit(i)
%          	2:      history{nite,1:6} = {A,B,fita,fitb,bestind,bestfit}
%   A0:             initial guess
%   goal:       If function value is below goal, iterations are stopped
%   nitemax:    maximum number of iterations allowed
%   mu:   	    Simulated annealing parameter (eg, 0.2, read below)
%
%   Call back functions to be provided by the user:
%   fitfun:	    fitness function
%   mutfun:	    mutation (change) function. Receives an individual and its
%               fitness and returns a modified individual
%   prifun:	    prints individual
%
%Outputs:
%   bestind:    best individual from the last generation
%   bestfit:    fitness value of best individual from the last population
%   nite:       number of iterations performed
%   lastpop:    population of last generation
%   lastfit:    fitness values of last population
%   history:    array with saved history array
%
%About mu parameter (from
%http://sourceforge.net/p/sbsi/discussion/1048774/thread/6caecc11/)
%High values of Mu mean that the probability of accepting a worse solution
%from one iteration to the next is low, and hence the solution is likely 
%to head straight to a local minimum. Conversely, low values of Mu mean 
%that the probability of accepting a worse solution from one iteration to 
%the next is high, and hence the solution is able to jump out of local 
%minima but may also not be able to converge at all. So finding a good 
%value of Mu is a tradeoff between these extremes and the absolute value 
%needs to be empirically determined for different models.

% Get options
if isfield(opts,'ninfo'), ninfo = opts.ninfo; else, ninfo = 1; end
if isfield(opts,'label'), label = opts.label; else, label = 0; end
if isfield(opts,'nhist'), nhist = opts.nhist; else, nhist = 1; end

% Create history array
history = [];

% Preprocessing
A = A0; % Set initial guess
fitA = fitfun(A); % Compute initial guess fitness
bestind = A; % Best individual
bestfit = fitA; % Best fitness

% Iterate until convergence or max iterations
for ite=1:nitemax
    
    % Mutate A into a new individual B
    B = mutfun(A,fitA);
    
    % Compute B fitness
    fitB = fitfun(B);

    % Save best individual/fitness
    if fitB<bestfit % New individual B has better fitness than old A
        bestind = B; % B is now best individual
        bestfit = fitB; % B has now best fitness
    end
    
    % Save history
    if nhist>1 % Save full history {A,B,fita,fitb}
        history{ite,1} = A; %#ok
        history{ite,2} = B; %#ok
        history{ite,3} = fitA; %#ok
        history{ite,4} = fitB; %#ok
        history{ite,5} = bestind; %#ok
        history{ite,6} = bestfit; %#ok
    elseif nhist>0 % Save best fitness only
        history(ite) = bestfit; %#ok
    end

    % Check if reached target fitness or max iterations
    if bestfit<goal || ite>=nitemax % Target achieved
        
        % Save last iteration data
        nite = ite; % Save current generation index
        lastpop = {A,B}; % Save last population
        lastfit = [fitA,fitB]; % Save last fitness values
        
        % Show info
        if ninfo>0
            fprintf('ASA label=%d nite=%2d fitbest=%f',label,nite,bestfit);
            if ~isempty(prifun), fprintf(' best='); prifun(bestind); end
            if bestfit<goal % Goal achieved
                fprintf(' goal=%e achieved, leaving\n',goal);
            else % Maximum generations reached (goal not achieved)
                fprintf(' max. iterations reached, leaving\n');
            end
        end
        
        % Stop iterating
        break;
        
    end
    
    % Print extended info
    if ninfo>1
        fprintf('ASA label=%d nite=%2d fitbest=%f',label,ite,bestfit);
        if ~isempty(prifun), fprintf(' best='); prifun(bestind); end
        if ninfo<=2, fprintf('\n'); end
    end

    % Compute fitness difference between B and A
    deltafit = fitB - fitA;
    
    % Simulate thermal transition probability (~ exp(-Delta/(k_b*T))
    probability = exp(-deltafit / (mu * abs(fitA)));
    
    % Print extra info
    if ninfo>2
        fprintf([' fitB=%8.2e deltafit=%+8.2e ', ...
            'deltafit/fit=%+8.2e probability=%+8.2e'], ...
            fitB,deltafit,deltafit/abs(fitA),probability);
    end
    
    % Jump from A to B (randomly with thermal transition probability)
    if rand<probability % Jump
        A = B; % Jump from A to B
        fitA = fitB; % Update fitness value
        if ninfo>2 % Print extra info
            fprintf(' jump '); 
            if deltafit<=0, fprintf('>= \n'); % Fitness improved
            else, fprintf('< \n'); % Fitness worsened
            end
        end
    else % Do not jump
        if ninfo>2 % Print extra info
            fprintf(' = \n'); % Fitness does not change
        end
    end
    
end

end

