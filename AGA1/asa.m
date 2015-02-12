function [ vA,fita ] = asa( npr,label, A0, nitemax, mu, goal, ffit, fnei,prifun )
% Iterates to find mimumum of a function using Simulated Annealing
% (c) Manel Soria - ETSEIAT
% npr:      iteration control; prints every ninfo iterations
% label:    integer number that precedes the prints in case output is to be
%           filtered    % A0: initial guess
% nitemax:  maxim number of iteracions
% mu:       Simulated annealing parameter (eg, 0.2, read below)
% goal:     If function value is below goal, iterations are stopped
%
% Call back functions to be provided by user:    
% fit       fitness function
% fnei      mutation (change) function. Receives an individual and its fitness
%           and returns a changed individual
% prifun:   Prints individual
%
% Returns
% vA        Best solution
% fita      Fitness of the best solution
%
% About mu parameter (from
% http://sourceforge.net/p/sbsi/discussion/1048774/thread/6caecc11/)
% High values of Mu mean that the probability of accepting a worse solution
% from one iteration to the next  is low, and hence the solution is likely 
% to head straight to a local minimum. Conversely,  low values of Mu mean 
% that the probability of accepting a worse solution from one iteration to 
% the next  is high, and hence the solution is able to jump out of local 
% minima but may also not be able to converge at all. So finding a good 
% value of Mu is a tradeoff between these extremes and the absolute value 
% needs to be empirically determined for different models. 

    info=0;

    A=A0;

    fita=ffit(A);
    best=A;
    fitbest=fita;    
    nite=1;
    while fita>goal && nite<=nitemax
        ppp=nite==1 || mod(nite,npr)==0;
        
        if ppp
            fprintf('SA label=%d nite=%d ',label,nite);
            prifun(A);
            fprintf(' fitA=%f fitbest=%f \n',fita,fitbest);            
        end
        
        if fita<goal
            break;
        end
        B=fnei(A,fita);        
        fitb=ffit(B);
        
        if fitb<fitbest % keep a copy of the best (just in case...)
            best=B;
            fitbest=fitb;
        end
        deltafit=fitb-fita;
        probability=exp(-deltafit/(abs(fita)*mu));
        if ppp && info
            fprintf('   fitb=%8.2e deltafit=%+8.2e deltafit/fit=%+8.2e probability=%+8.2e  ',...
                        fitb,      deltafit,       deltafit/abs(fita),probability);  
        end      
        if rand<probability 
            A=B;
            fita=fitb;
            if ppp && info
                fprintf('jump '); 
                if deltafit<=0
                    fprintf('>= \n'); 
                else
                    fprintf('< \n'); 
                end
            end
        else
            if ppp && info, fprintf('= \n'); end
        end
        nite=nite+1;
    end
    if fitbest<fita
        A=best;
        fita=fitbest;
    end
    vA=A;
end

