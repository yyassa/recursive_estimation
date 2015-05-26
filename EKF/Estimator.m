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
    
    % Initialize position/orientation with zero & wheel radius with W_0
    posEst = [0 0];
    oriEst = 0;
    radiusEst = knownConst.NominalWheelRadius;
    
    % Initialize variances with uniform dist variance.
    posVar = 1/12*(2*knownConst.TranslationStartBound)^2*[1 1];
    oriVar = 1/12*(2*knownConst.RotationStartBound)^2;
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
 
tspan = [estState.last_tm, tm];

% Define: xP - [20x1] state/covariance vector: 
% xP(1:4) ~ states 
% xP(5:20) ~ stacked elements of P

% Kinematic equations as given in problem
s_v = @(xP) xP(4)*actuate(1);
s_t = @(xP) s_v(xP)*cos(actuate(2));
s_r = @(xP) -1/knownConst.WheelBase * s_v(xP) * sin(actuate(2));

% Differentiate system dynamics wrt states
A = @(t,xP) [0 0 -s_t(xP)*sin(xP(3)) actuate(1)*cos(actuate(2))*cos(xP(3)); ...
             0 0  s_t(xP)*cos(xP(3)) actuate(1)*cos(actuate(2))*sin(xP(3)); ...
             0 0                    0            -1/knownConst.WheelBase*actuate(1)*sin(actuate(2)); ...
             0 0                    0                                          0];

if(designPart == 1)
    
    L = eye(4);
    Q = diag([0.1,0.1,0.001,0]);
    
    % Stacked ODE
    q = @(t,xP) [ % state dynamics
                   s_t(xP)*cos(xP(3)); ...
                   s_t(xP)*sin(xP(3)); ...
                   s_r(xP); ...
                   0; ...
                   % covariance dynamic
                   reshape( ... %vectorize 4x4 covariance to 16x1
                        A(t,xP)*reshape(xP(5:20),[4 4]) ... 
                        + reshape(xP(5:20),[4 4])*(A(t,xP))' + L*Q*L', ...
                    [16 1])];
    
    % Initialize and solve ODE
    states0 = estState.states(1:4);
    P0 = reshape(estState.P, [16 1]);
    [~, sol] = ode45(q,tspan,[states0;P0]);

    x_p = sol(end,1);
    y_p = sol(end,2);
    r_p = sol(end,3);
    W_p = sol(end,4);
    states_p = [x_p; y_p; r_p; W_p];
    P_p = reshape(sol(end,5:20),[4 4]);

elseif (designPart == 2)
    % With the new process noise model we do have only 2 noise components
    % (v_v and v_r) instead of 4, but with known variances Q_v and Q_r. The 
    % linearization of the dynamics wrt these noiseterms leads to:
    %   modified L(t) matrix with dimension [4x2],
    %   modified Q matrix with dimension [2x2]
    
    q = @(t,x) [s_t(x)*cos(x(3)); s_t(x)*sin(x(3)); s_r(x); 0];
    states0 = estState.states(1:4);
    [t_vector, sol] = ode45(q,tspan,states0);

    x_p = sol(end,1);
    y_p = sol(end,2);
    r_p = sol(end,3);
    W_p = sol(end,4);

states_p = [x_p; y_p; r_p; W_p];
    
    Q = diag([knownConst.VelocityInputPSD, knownConst.AngleInputPSD]); % dimension [2x2]
    
    qP = @(t,P) reshape(A(t,t_vector,sol,actuate, knownConst.WheelBase)*reshape(P,[4 4]) + reshape(P,[4 4])*(A(t,t_vector,sol,actuate, knownConst.WheelBase))' + L(t,t_vector,sol,actuate, knownConst.WheelBase)*Q*(L(t,t_vector,sol,actuate, knownConst.WheelBase))', [16 1]); 
    P0 = reshape(estState.P, [16 1]);
    [~,solP] = ode45(qP, tspan,P0);
    P_p = reshape(solP(end,:),[4 4]);
    
end



%% S2: Measurement update

% Incase of no measurement, set posterior estimate equal to prior update
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
    
    % Incase an orientation measurement follows, trick it into using the
    % distance measurement updated estimate as its prior.
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


posEst = [states_m(1) states_m(2)];
oriEst = states_m(3);
posVar = [P_m(1,1) P_m(2,2)];
oriVar = P_m(3,3);
radiusEst = states_m(4);
radiusVar = P_m(4,4);

estState.states = states_m;
estState.P = P_m;

% remember timestamp for next iteration for the cont. ode solver
estState.last_tm = tm;
end

% function A = A(t, t_vector, sol, actuate, B)
% % A has dimension [4x4]
% 
%     %find the closest timestamp to t from the first ode's time vector
%     %(state and variance propagation use seperate ode funcs -> different time vectors) 
%     [~,ind] = min(abs(t-t_vector));
% 
%     s_v = sol(ind,4)*actuate(1);
%     s_t = s_v*cos(actuate(2));
%     % differentiate system dynamics wrt states
%     A = [0 0 -s_t*sin(sol(ind,3)) actuate(1)*cos(actuate(2))*cos(sol(ind,3)); ...
%          0 0  s_t*cos(sol(ind,3)) actuate(1)*cos(actuate(2))*sin(sol(ind,3)); ...
%          0 0                    0            -1/B*actuate(1)*sin(actuate(2)); ...
%          0 0                    0                                          0];
% end

function L = L(t, t_vector, sol, actuate, B)
% L has dimension [4x2]
    
    %find the closest timestamp to t from the first ode's time vector
    %(state and variance propagation use seperate ode funcs -> different time vectors) 
    [~,ind] = min(abs(t-t_vector));
    
    % differentiate 
    L = [sol(ind,4)*actuate(1)*cos(actuate(2))*cos(sol(ind,3)) -sol(ind,4)*actuate(1)*sin(actuate(2))*cos(sol(ind,3)); ...
         sol(ind,4)*actuate(1)*cos(actuate(2))*sin(sol(ind,3)) -sol(ind,4)*actuate(1)*sin(actuate(2))*sin(sol(ind,3)); ...
         -(1/B)*sol(ind,4)*actuate(1)*sin(actuate(2))          -(1/B)*sol(ind,4)*actuate(1)*cos(actuate(2))          ; ...
         0                                                     0                                                    ];

end
