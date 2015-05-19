function [postParticles] = Estimator(prevPostParticles, sens, act, init)
% [postParticles] = Estimator(prevPostParticles, sens, act, init)
%
% The estimator function. The function will be called in two different
% modes: If init==1, the estimator is initialized. If init == 0, the
% estimator does an iteration for a single sample time interval Ts (KC.ts)
% using the previous posterior particles passed to the estimator in
% prevPostParticles and the sensor measurements and control inputs.
%
% You must edit this function.
%
% Inputs:
%   prevPostParticles   previous posterior particles at discrete time k-1,
%                       which corresponds to continuous time t = (k-1)*Ts
%                       The variable is a struct whose fields are arrays
%                       that correspond to the posterior particle states.
%                       The fields are: (N is number of particles)
%                       .x = (2xN) array with the x-locations (metres)
%                       .y = (2xN) array with the y-locations (metres)
%                       .h = (2xN) array with the headings (radians)
%                       The first row in the arrays corresponds to robot A.
%                       The second row corresponds to robot B.
%
%   sens                Sensor measurements at discrete time k (t = k*Ts),
%                       [4x1]-array, an Inf entry indicates no measurement
%                       of the corresponding sensor.
%                       sens(1): distance reported by sensor 1 (metres)
%                       sens(2): distance reported by sensor 2 (metres)
%                       sens(3): distance reported by sensor 3 (metres)
%                       sens(4): distance reported by sensor 4 (metres)
%
%   act                 Control inputs u at discrete time k-1, which are
%                       constant during a time interval Ts:
%                       u(t) = u(k-1) for (k-1)*Ts <= t < k*Ts
%                       [2x1]-array:
%                       act(1): velocity of robot A, u_A(k-1) (metres/second)
%                       act(2): velocity of robot B, u_B(k-1) (metres/second)
%
%   init                Boolean variable indicating wheter the estimator
%                       should be initialized (init = 1) or if a regular
%                       estimator update should be performed (init = 0).
%                       OPTIONAL ARGUMENT. By default, init = 0.
%
% Outputs:
%   postParticles       Posterior particles at discrete time k, which
%                       corresponds to the continuous time t = k*Ts.
%                       The variable is a struct whose fields are arrays
%                       that correspond to the posterior particle states.
%                       The fields are: (N is number of particles)
%                       .x = (2xN) array with the x-locations (metres)
%                       .y = (2xN) array with the y-locations (metres)
%                       .h = (2xN) array with the headings (radians)
%                       The first row in the arrays corresponds to robot A.
%                       The second row corresponds to robot B.
%
% Class:
% Recursive Estimation
% Spring 2015
% Programming Exercise 2
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Mark Mueller
% mwm@ethz.ch

% Check if init argument was passed to estimator:
if(nargin < 4)
    % if not, set to default value:
    init = 0;
end

%% Mode 1: Initialization
% Set number of particles:
N = 1000; % obviously, you will need more particles than 10.
if (init)
    % Do the initialization of your estimator here!
    % These particles are the posterior particles at discrete time k = 0
    % which will be fed into your estimator again at k = 1
    % Replace the following:
    postParticles.x = [KC.L*ones(1,N);zeros(1,N)];
    postParticles.y = (randi(2,2,N)-1).*KC.L;
    B = postParticles.y/KC.L;
    postParticles.h = [ (pi/2).*rand(1,N)+pi/2-B(1,:)*3*pi/2; ...
                        (pi/2).*rand(1,N)-B(2,:)*pi/2];
    return;
end % end init

%% Mode 2: Estimator iteration.
% If init = 0, we perform a regular update of the estimator.

% S1:Process update
    priorParticles.x = prevPostParticles.x + diag(act)*cos(prevPostParticles.h)*KC.ts;
    priorParticles.y = prevPostParticles.y + diag(act)*sin(prevPostParticles.h)*KC.ts;
    priorParticles.h = prevPostParticles.h;
    
    %Lower wall B1
    B1 = priorParticles.y <= 0 & diag(act)*sin(priorParticles.h)<0;
    priorParticles.y(B1) = 0;
    v = (rand(2,N)).^(1/3)*KC.vbar.*(randi(2,2,N)*2-3);
    alpha = -priorParticles.h.*(1+v); 
    priorParticles.h(B1) = alpha(B1);
    
    %Lower wall B1
    B1 = priorParticles.y <= 0 & diag(act)*sin(priorParticles.h)<0;
    priorParticles.y(B1) = 0;
    v = (rand(2,N)).^(1/3)*KC.vbar.*(randi(2,2,N)*2-3);
    alpha = -priorParticles.h.*(1+v); 
    priorParticles.h(B1) = alpha(B1);
    
    %Right wall B2
    B2 = priorParticles.x >= KC.L & diag(act)*cos(priorParticles.h)>0;
    priorParticles.x(B2) = KC.L;
    v = (rand(2,N)).^(1/3)*KC.vbar.*(randi(2,2,N)*2-3);
    alpha = (pi/2 - priorParticles.h).*(1+v); 
    priorParticles.h(B2) = pi/2 + alpha(B2);
    
    %Upper wall B3
    B3 = priorParticles.y >= KC.L & diag(act)*sin(priorParticles.h)>0;
    priorParticles.y(B3) = KC.L;
    v = (rand(2,N)).^(1/3)*KC.vbar.*(randi(2,2,N)*2-3);
    alpha = priorParticles.h.*(1+v); 
    priorParticles.h(B3) = -alpha(B3);
    
    %Left wall B4
    B4 = priorParticles.x <= 0 & diag(act)*cos(priorParticles.h)<0;
    priorParticles.x(B4) = 0;
    v = (rand(2,N)).^(1/3)*KC.vbar.*(randi(2,2,N)*2-3);
    alpha = (priorParticles.h-pi/2).*(1+v); 
    priorParticles.h(B4) = pi/2-alpha(B4);
    
    priorParticles.h = mod(priorParticles.h+pi,2*pi)-pi;
    
% S2:Measurement update
    postParticles = priorParticles;
    
    if(isfinite(sens(1)))
         d = sqrt((priorParticles.x - KC.L).^2 + priorParticles.y.^2);
         %lik_* has to be looped through all particles.
         lik_cor = makedist('Triangular','a',sens(1) - d(1,:) - KC.wbar,'b',sens(1) - d(1,:),'c',sens(1) - d(1,:) + KC.wbar);
         lik_false = makedist('Triangular','a',sens(1) - d(2,:) - KC.wbar,'b',sens(1) - d(2,:),'c',sens(1) - d(1,:) + KC.wbar);
         beta = sbar*random(lik_false,1,N)+(1-sbar)*random(lik_cor,1,N);
         alpha = 1/sum(beta);
         beta = alpha*beta;
         beta_kum = cumsum(beta);
       
         for j = 1:N;
            r = rand();
            postParticles.x = 0;
         end
    end
    
% Replace the following:
% postParticles.x = zeros(2,N);
% postParticles.y = zeros(2,N);
% postParticles.h = zeros(2,N);

end % end estimator

