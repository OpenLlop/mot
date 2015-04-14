function [ bestind, bestfit, nite, history ] = asa ( opts, ...
    A0, nitemax, mu, goal, ...
    fitfun, mutfun, prifun )
% Iterates to find mimumum of a function using Simulated Annealing
%
%Programmers:   Manel Soria         (UPC/ETSEIAT)
%               David de la Torre	(UPC/ETSEIAT)
%Date:          14/04/2015
%Revision:      2
%
%Usage:         [bestind, bestfit, nite, history] = asa ( opts, ...
%                   A0, nitemax, mu, goal, ...
%                   fitfun, mutfun, prifun )
%
%Inputs:
%   opts:       function control parameters [struct]
%       ninfo:  iteration control; prints every ninfo iterations
%       einfo:  prints extended information every ninfo iterations [1,0]
%       label:  integer number that precedes the prints in case output is
%               to be filtered
%       fhist:  saved history level (0=none, 1=just fitness, 2=all data)
%                   0: history = []
%                   1: history(nite) = bestfit(i)
%                   2: history{nite,1:6} = {A,B,fita,fitb,bestind,bestfit}
%   A0:         initial guess
%   nitemax:    maxim number of iteracions
%   mu:         Simulated annealing parameter (eg, 0.2, read below)
%   goal:       If function value is below goal, iterations are stopped
%
%   Call back functions to be provided by user:
%   fitfun:     fitness function
%   mutfun:     mutation (change) function. Receives an individual and its
%               fitness and returns a modified individual
%   prifun:     Prints individual
%
%Outputs:
%   bestind:    Best individual
%   bestfir:    Fitness of the best individual
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

% Set default options
if isfield(opts,'ninfo'), ninfo = opts.ninfo; else ninfo = 1; end;
if isfield(opts,'einfo'), einfo = opts.einfo; else einfo = 0; end;
if isfield(opts,'label'), label = opts.label; else label = 0; end;
if isfield(opts,'nhist'), nhist = opts.nhist; else nhist = 1; end;

% Declare history array, if required
if nhist>0, history = []; end;

% Preprocessing
A = A0; % Set initial guess
fitA = fitfun(A); % Compute initial guess fitness
bestind = A; % Best individual
bestfit = fitA; % Best fitness

% Iterate until convergence or max iterations
for nite=1:nitemax
    
    % Must print something on this iteration
    mustprint = nite==1 || mod(nite,ninfo)==0;
    
    % Print info
    if mustprint
        fprintf('SA label=%d nite=%d ',label,nite); % Print info
        prifun(A); % Print individual
        fprintf(' fitA=%f fitbest=%f \n',fitA,bestfit); % Print fitness
    end;

    % Check if reached target fitness or max iterations 
    if fitA<goal || nite>=nitemax % Target achieved; end simulation
        return; % Return from function
    end;
    
    % Mutate A into a new individual B
    B = mutfun(A,fitA);
    
    % Compute B fitness
    fitB = fitfun(B);

    % Save best individual/fitness
    if fitB<bestfit % New individual B has better fitness than old A
        bestind = B; % B is now best individual
        bestfit = fitB; % B has now best fitness
    end;
    
    % Save history
    if nhist>1 % Save full history {A,B,fita,fitb}
        history{nite,1} = A; %#ok
        history{nite,2} = B; %#ok
        history{nite,3} = fitA; %#ok
        history{nite,4} = fitB; %#ok
        history{nite,5} = bestind; %#ok
        history{nite,6} = bestfit; %#ok
    elseif nhist>0 % Save best fitness only
        history(nite) = bestfit; %#ok
    end;

    % Compute fitness difference between B and A
    deltafit = fitB - fitA;
    
    % Simulate thermal transition probability (~ exp(-Delta/(k_b*T))
    probability = exp(-deltafit / (mu * abs(fitA)));
    
    % Print extra info
    if mustprint && einfo
        fprintf(['   fitB=%8.2e deltafit=%+8.2e ', ...
            'deltafit/fit=%+8.2e probability=%+8.2e  '], ...
            fitB,deltafit,deltafit/abs(fitA),probability);
    end;
    
    % Jump from A to B (randomly with thermal transition probability)
    if rand<probability % Jump
        A = B; % Jump from A to B
        fitA = fitB; % Update fitness value
        if mustprint && einfo % Print extra info
            fprintf('jump '); 
            if deltafit<=0, fprintf('>= \n'); % Fitness improved
            else fprintf('< \n'); % Fitness worsened
            end;
        end;
    else % Do not jump
        if mustprint && einfo % Print extra info
            fprintf('= \n'); % Fitness does not change
        end;
    end;
    
end;

end

