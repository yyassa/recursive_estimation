function [posEst,oriEst,radiusEst, posVar,oriVar,radiusVar,estState] = Estimator(estState,actuate,sense,tm,knownConst,designPart)
% [posEst,oriEst,posVar,oriVar,baseEst,baseVar,estState] =
% 	Estimator(estState,actuate,sense,tm,knownConst,designPart)
%
% The estimator.
%
% The Estimator function shall be used for both estimator design parts; the
% input argument designPart is used to distinguish the two:
%   designPart==1  -> Part 1
%   designPart==2  -> Part 2
%
% The function will be called in two different modes:
% If tm==0, the estimator is initialized; otherwise the estimator does an
% iteration step (compute estimates for the time step k).
%
% Inputs:
%   estState        previous estimator state (time step k-1)
%                   May be defined by the user (for example as a struct).
%   actuate         control input u(k-1), [1x2]-vector
%                   actuate(1): u_v(k-1), drive wheel angular velocity
%                   actuate(2): u_r(k-1), drive wheel angle
%   sense           sensor measurements z(k), [1x2]-vector, INF if no
%                   measurement
%                   sense(1): z_d(k), distance measurement
%                   sense(2): z_r(k), orientation measurement
%   tm              time, scalar
%                   If tm==0 initialization, otherwise estimator
%                   iteration step.
%   knownConst      known constants (from KnownConstants.m)
%   designPart      variable to distinguish the estimator design part
%                       designPart==1  -> Part 1
%                       designPart==2  -> Part 2
%
% Outputs:
%   posEst          position estimate (time step k), [1x2]-vector
%                   posEst(1): x position estimate
%                   posEst(2): y position estimate
%   oriEst          orientation estimate (time step k), scalar
%   radiusEst       estimate of wheel radius W (time step k), scalar
%   posVar          variance of position estimate (time step k), [1x2]-vector
%                   posVar(1): x position variance
%                   posVar(2): y position variance
%   oriVar          variance of orientation estimate (time step k), scalar
%   radiusVar       variance of wheel radius estimate (time step k), scalar
%   estState        current estimator state (time step k)
%                   Will be input to this function at the next call.
%
%
% Class:
% Recursive Estimation
% Spring 2015
% Programming Exercise 1
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Michael Muehlebach
% michaemu@ethz.ch
%
% --
% Revision history
% [19.04.11, ST]    first version by Sebastian Trimpe
% [30.04.12, PR]    adapted version for spring 2012, added unknown wheel
%                   radius
% [06.05.13, MH]    2013 version
% [24.04.15, MM]    2015 version


%% Mode 1: Initialization
if (tm == 0)
    % Do the initialization of your estimator here!
    
    % Replace the following:
    posEst = [0 0];
    oriEst = 0;
    posVar = 1/12*(2*knownConst.TranslationStartBound)^2*[1 1];
    oriVar = 1/12*(2*knownConst.RotationStartBound)^2;
    radiusEst = knownConst.NominalWheelRadius;
    radiusVar = 1/12*(2*knownConst.WheelRadiusError)^2;  
    
    estState.states = [posEst(1); posEst(2); oriEst; radiusEst];
    estState.P = diag([posVar(1) posVar(2) oriVar radiusVar]);
    estState.last_tm = tm;
    return;
end


%% Mode 2: Estimator iteration.
% If we get this far tm is not equal to zero, and we are no longer
% initializing.  Run the estimator.

%% S1: Prior update
 
% propagate States x_p
tspan = [estState.last_tm, tm];

s_v = @(x) x(4)*actuate(1);
s_t = @(x) s_v(x)*cos(actuate(2));
s_r = @(x) -1/knownConst.WheelBase * s_v(x) * sin(actuate(2));

q = @(t,x) [s_t(x)*cos(x(3)); s_t(x)*sin(x(3)); s_r(x); 0];
x0 = estState.states(1:4);
[t_vector, sol] = ode45(q,tspan,x0);

x_p = sol(end,1);
y_p = sol(end,2);
r_p = sol(end,3);
W_p = sol(end,4);

states_p = [x_p; y_p; r_p; W_p];


% propagate Covariances P_p
if(designPart == 1)
    L_1 = eye(4);
    Q = diag([0.1,0.1,0.001,0]);

    qP = @(t,P) reshape(A(t,t_vector,sol,actuate, knownConst.WheelBase)*reshape(P,[4 4]) + reshape(P,[4 4])*(A(t,t_vector,sol,actuate, knownConst.WheelBase))' + L_1*Q*L_1', [16 1]); 
    P0 = reshape(estState.P, [16 1]);
    [~,solP] = ode45(qP, tspan,P0); 

elseif (designPart == 2)

    % Because we have two seperate noises now that cannot be simply added to
    % eachother, we need to introduce one more L and Q and adapt the
    % derivation. 
    %
    % For better readability, we can make them three seperate Ls and Qs (for
    % Q_v,Q_r and the Wheel Bias)..using 2 or 3 matrices here is mathematically
    % equivalent, as we mostly fill the last one with zeros!
    %
    % (!) - Unfortunately, the Ls here are time dependent, so it gets shitty
    %       again with plugging it into the ODE.

%     %L_v: differentiated w.r.t v_v (see function)
%     Q_v = eye(4).*knownConst.AngleInputPSD;
%     Q_v(4,4) = 0;
% 
%     %L_r: differentiated w.r.t v_r (see function)
%     Q_r = eye(4).*knownConst.VelocityInputPSD;
%     Q_r(4,4) = 0;
% 
%     L_W = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 1];
%     Q_W = diag([0,0,0,0.0001]);
% 
%     % qP = APA' + LQL' + LQL' + LQL' (jeweils L_v,L_r,L_W und Q_v,Q_r,Q_W)
%     qP = @(t,P) reshape(A(t,t_vector,sol,actuate, knownConst.WheelBase)*reshape(P,[4 4]) + reshape(P,[4 4])*(A(t,t_vector,sol,actuate, knownConst.WheelBase))' + L_v(t,t_vector,sol,actuate, knownConst.WheelBase)*Q_v*(L_v(t,t_vector,sol,actuate, knownConst.WheelBase))' + L_r(t,t_vector,sol,actuate, knownConst.WheelBase)*Q_r*(L_r(t,t_vector,sol,actuate, knownConst.WheelBase))' + L_W*Q_W*L_W', [16 1]); 
%     P0 = reshape(estState.P, [16 1]);
%     [~,solP] = ode45(qP, tspan,P0); 
    
    Q = diag([knownConst.VelocityInputPSD, knownConst.AngleInputPSD, 0]); % dimension [3x3]
    
    qP = @(t,P) reshape(A(t,t_vector,sol,actuate, knownConst.WheelBase)*reshape(P,[4 4]) + reshape(P,[4 4])*(A(t,t_vector,sol,actuate, knownConst.WheelBase))' + L(t,t_vector,sol,actuate, knownConst.WheelBase)*Q*(L(t,t_vector,sol,actuate, knownConst.WheelBase))', [16 1]); 
    P0 = reshape(estState.P, [16 1]);
    [~,solP] = ode45(qP, tspan,P0);

end

P_p = reshape(solP(end,:),[4 4]);

%% S2: Measurement update
states_m = states_p;
P_m = P_p;

% Distance Measurement
if(isfinite(sense(1)))

    H = 1/norm([x_p,y_p])*[ x_p y_p 0 0];
    M = 1;
    R = 1/6*(knownConst.DistNoise)^2;
    
    K = P_p*H'*inv( H*P_p*H' + M*R*M');
    states_m = states_p + K*(sense(1) - norm([x_p,y_p]));
    P_m = (eye(4) - K*H)*P_p;
    
    states_p = states_m;
    P_p = P_m;
end
% Orientation Measurement
if(isfinite(sense(2)))
   
    H = [ 0 0 1 0];
    M = 1;
    R = knownConst.CompassNoise;
    
    K = P_p*H'*inv( H*P_p*H' + M*R*M');
    states_m = states_p + K*(sense(2) - states_p(3));
    P_m = (eye(4) - K*H)*P_p;
end


% Replace the following:
posEst = [states_m(1) states_m(2)];
oriEst = states_m(3);
posVar = [P_m(1,1) P_m(2,2)];
oriVar = P_m(3,3);
radiusEst = states_m(4);
radiusVar = P_m(4,4);

estState.states = [posEst(1); posEst(2); oriEst; radiusEst];
estState.P = diag([posVar(1) posVar(2) oriVar radiusVar]);
estState.last_tm = tm;
end

function A = A(t, t_vector, sol, actuate, B)
% A has dimension [4x4]
    [~,ind] = min(abs(t-t_vector));

    s_v = sol(ind,4)*actuate(1);
    s_t = s_v*cos(actuate(2));   
    A = [0 0 -s_t*sin(sol(ind,3)) actuate(1)*cos(actuate(2))*cos(sol(ind,3)); ...
         0 0 s_t*cos(sol(ind,3)) actuate(1)*cos(actuate(2))*sin(sol(ind,3)); ...
         0 0 0 -1/B*actuate(1)*sin(actuate(2)); ...
         0 0 0 0];
end

function L = L(t, t_vector, sol, actuate, B)
% L has dimension [4x3]
% --> Matrix Q needs to be of size [3x3]
    [~,ind] = min(abs(t-t_vector));
    
    L = [sol(ind,4)*actuate(1)*cos(actuate(2))*cos(sol(ind,3)) -sol(ind,4)*actuate(1)*sin(actuate(2))*cos(sol(ind,3)) 0; ...
         sol(ind,4)*actuate(1)*cos(actuate(2))*sin(sol(ind,3)) -sol(ind,4)*actuate(1)*sin(actuate(2))*sin(sol(ind,3)) 0; ...
         -(1/B)*sol(ind,4)*actuate(1)*sin(actuate(2))          -(1/B)*sol(ind,4)*actuate(1)*cos(actuate(2))           0; ...
         0                                                     0                                                      1];

end

function L = L_v(t,t_vector,sol,actuate,B)
    [~,ind] = min(abs(t-t_vector));
    
    L = diag([sol(ind,4)*actuate(1)*cos(actuate(2))*cos(sol(ind,3)) ...
              sol(ind,4)*actuate(1)*cos(actuate(2))*sin(sol(ind,3)) ...
              -(1/B)*sol(ind,4)*actuate(1)*sin(actuate(2))          ...
              0]);
    % probably instead of sol(ind,4) here we could use the constant last
    % estimate of zero because i doesn't change over time..just being
    % super-consistent.
end

function L = L_r(t,t_vector,sol,actuate,B)
    [~,ind] = min(abs(t-t_vector));
    
    L = diag([-sol(ind,4)*actuate(1)*sin(actuate(2))*cos(sol(ind,3)) ...
              -sol(ind,4)*actuate(1)*sin(actuate(2))*sin(sol(ind,3)) ...
              -(1/B)*sol(ind,4)*actuate(1)*cos(actuate(2))          ...
              0]);
    
    % probably instead of sol(ind,4) here we could use the constant last
    % estimate of zero because i doesn't change over time..just being
    % super-consistent.
end
