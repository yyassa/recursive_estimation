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
N = 2000;

if (init)
    % These particles are the posterior particles at discrete time k = 0
    % which will be fed into your estimator again at k = 1
    postParticles.x = [KC.L*ones(1,N);zeros(1,N)];
    postParticles.y = (randi(2,2,N)-1).*KC.L;
    B = postParticles.y/KC.L;
    postParticles.h = [ (pi/2).*rand(1,N)+pi/2-B(1,:)*3*pi/2; ...
                        (pi/2).*rand(1,N)-B(2,:)*pi/2];
    return;
end % end init

%% Mode 2: Estimator iteration.
% If init = 0, we perform a regular update of the estimator.

% ------------------------ S1: Process update --------------------------- %

% Propagate particles without walls:
    priorParticles.x = prevPostParticles.x + diag(act)*cos(prevPostParticles.h)*KC.ts;
    priorParticles.y = prevPostParticles.y + diag(act)*sin(prevPostParticles.h)*KC.ts;
    priorParticles.h = prevPostParticles.h;
    
% Lower wall B1
    v = (rand(2,N)).^(1/3)*KC.vbar.*(randi(2,2,N)*2-3);
    % case a: h >= -90deg
    B1a = priorParticles.y <= 0 & diag(act)*sin(priorParticles.h)<0 & priorParticles.h >= -pi/2;
    priorParticles.y(B1a) = 0;
    alpha = -priorParticles.h.*(1+v); 
    priorParticles.h(B1a) = alpha(B1a);
    % case b: h < -90deg
    B1b = priorParticles.y <= 0 & diag(act)*sin(priorParticles.h)<0 & priorParticles.h < -pi/2;
    priorParticles.y(B1b) = 0;
    alpha = (pi + priorParticles.h).*(1+v); 
    priorParticles.h(B1b) = pi - alpha(B1b);
    
% Right wall B2
    v = (rand(2,N)).^(1/3)*KC.vbar.*(randi(2,2,N)*2-3);
    % case a: h >= 0
    B2a = priorParticles.x >= KC.L & diag(act)*cos(priorParticles.h)>0 & priorParticles.h >= 0;
    priorParticles.x(B2a) = KC.L;
    alpha = (pi/2 - priorParticles.h).*(1+v); 
    priorParticles.h(B2a) = pi/2 + alpha(B2a);
    % case b: h < 0
    B2b = priorParticles.x >= KC.L & diag(act)*cos(priorParticles.h)>0 & priorParticles.h < 0;
    priorParticles.x(B2b) = KC.L;
    alpha = (pi/2 + priorParticles.h).*(1+v); 
    priorParticles.h(B2b) = -pi/2 - alpha(B2b);
    
% Upper wall B3
    v = (rand(2,N)).^(1/3)*KC.vbar.*(randi(2,2,N)*2-3);
    % case a: h <= 90deg
    B3a = priorParticles.y >= KC.L & diag(act)*sin(priorParticles.h)>0 & priorParticles.h <= pi/2;
    priorParticles.y(B3a) = KC.L;
    alpha = priorParticles.h.*(1+v); 
    priorParticles.h(B3a) = -alpha(B3a);
    % case a: h > 90deg
    B3b = priorParticles.y >= KC.L & diag(act)*sin(priorParticles.h)>0 & priorParticles.h > pi/2;
    priorParticles.y(B3b) = KC.L;
    alpha = (pi - priorParticles.h).*(1+v); 
    priorParticles.h(B3b) = -pi + alpha(B3b);
    
% Left wall B4
    v = (rand(2,N)).^(1/3)*KC.vbar.*(randi(2,2,N)*2-3);
    % case a: h <= 0
    B4a = priorParticles.x <= 0 & diag(act)*cos(priorParticles.h)<0 & priorParticles.h <= 0;
    priorParticles.x(B4a) = 0;
    alpha = (pi/2 - priorParticles.h).*(1+v); 
    priorParticles.h(B4a) = -pi/2 + alpha(B4a);
    % case a: h > 0
    B4b = priorParticles.x <= 0 & diag(act)*cos(priorParticles.h)<0 & priorParticles.h > 0;
    priorParticles.x(B4b) = 0;
    alpha = (priorParticles.h-pi/2).*(1+v); 
    priorParticles.h(B4b) = pi/2 - alpha(B4b);
    
% Make sure the heading is between -pi and pi
    priorParticles.h = mod(priorParticles.h+pi,2*pi)-pi;
    
% ---------------------- S2: Measurement update ------------------------- %
% If there is no measurement use Prior:
postParticles = priorParticles;

% Save Sensorlocations:
x_sens = [KC.L KC.L 0    0]; % position of sensors in x
y_sens = [0    KC.L KC.L 0]; % position of sensors in y

for i = 1:4    % for all Sensors
    if(isfinite(sens(i)))  % if a measurement is published do
        % preallocate
        beta = zeros(1,N);
        
        for k = 1:N % for all Particles do
            
            % calculate distances to both robots
            d_1A = sqrt((priorParticles.x(1,k) - x_sens(i)).^2 + (priorParticles.y(1,k) - y_sens(i)).^2);
            d_1B = sqrt((priorParticles.x(2,k) - x_sens(i)).^2 + (priorParticles.y(2,k) - y_sens(i)).^2);
            
            if i<3 % Sensors 1 and 2
                
                % probability according to triangular distribution (linear
                % inside region [-wbar,wbar], 0 otherwise)
                lik_correct = 0;
                lik_false   = 0;
                if abs(sens(i)-d_1A) < KC.wbar
                    lik_correct = 1/KC.wbar - abs(sens(i)-d_1A)/KC.wbar^2;
                end
                if abs(sens(i)-d_1B) < KC.wbar
                    lik_false = 1/KC.wbar - abs(sens(i)-d_1B)/KC.wbar^2;
                end
                
            else % Sensors 3 and 4
                
                % probability according to triangular distribution (linear
                % inside region [-wbar,wbar], 0 otherwise)
                lik_correct = 0;
                lik_false   = 0;
                if abs(sens(i)-d_1B) < KC.wbar
                    lik_correct = 1/KC.wbar - abs(sens(i)-d_1B)/KC.wbar^2;
                end
                if abs(sens(i)-d_1A) < KC.wbar
                    lik_false = 1/KC.wbar - abs(sens(i)-d_1A)/KC.wbar^2;
                end
            end
            
            % calculate likelihood according to total probability theorem
            beta(k) = KC.sbar*lik_false + (1-KC.sbar)*lik_correct;
        end
        
        % normalise beta
        alpha = 1/sum(beta);
        beta = alpha*beta;
        
        % cummulate (beta_kum contains values from 0 to 1 in increasing order)
        beta_kum = cumsum(beta);
       
        % Reinizialise randomly if filter runs into numerical issues:
        if any(isnan(beta_kum))
            postParticles.x = KC.L*rand(2,N);
            postParticles.y = KC.L*rand(2,N);
            postParticles.h = 2*pi*rand(2,N)-pi;
            return;
        end
        
        % Resampling:
        for j = 1:N;
           r = rand();
           nbar = find(beta_kum >= r,1);
           postParticles.x(:,j) = priorParticles.x(:,nbar);
           postParticles.y(:,j) = priorParticles.y(:,nbar);
           postParticles.h(:,j) = priorParticles.h(:,nbar);
        end
        
        % If there is more than one measurement at the same time:
        priorParticles = postParticles;
    end
end
    
% Roughening

K = 0.01;
E_xy =sqrt(2)*KC.L;    % maximum inter-sample variability
E_h = pi;              % maximum inter-sample variability

sigma_xy = K*E_xy*N^(1/6);
sigma_h = K*E_h*N^(1/6);


delta_x = randn(2,N)*sigma_xy^2;
delta_y = randn(2,N)*sigma_xy^2;
delta_h = randn(2,N)*sigma_h^2;

postParticles.x = limit(postParticles.x + delta_x,0,KC.L);
postParticles.y = limit(postParticles.y + delta_y,0,KC.L);
postParticles.h = mod(postParticles.h + delta_h + pi,2*pi) - pi; %limit(postParticles.h + delta_h,-pi,pi); 

end % end estimator

function a = limit(a, lowerbound, upperbound)
    a(a < lowerbound) = lowerbound;
    a(a > upperbound) = upperbound;
end

