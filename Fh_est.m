function [F_h, F_h_rel, Imp_h, F_h_adj, F_h_adj_rel, Imp_h_adj] = Fh_est(t_c, t_f, v, m, L)
%% Simple hGRF Estimation

% Author:
% Geoffrey Burns
% University of Michigan

% Inputs:
% t_f -- flight time (s)
% t_c -- contact time (s)
% v   -- subject's average running speed (m/s)
% m   -- subject's mass (kg)
% L   -- subject's leg length (m)
%        Approximate leg length via height: L = 0.53*h
%        where h is the subject's height (m) per Winter (2005)

% Outputs:
% F_h         -- peak hGRF (N) (negative in braking/positive in propulsion)
% F_h_rel     -- peak relative hGRF (BW) (negative in braking/positive in propulsion)
% Imp_h       -- hGRF impulse (N-s) (negative in braking/positive in propulsion)
% F_h_adj     -- speed-corrected peak hGRF
% F_h_adj_rel -- speed-corrected relative peak hGRF (BW)
% Imp_h_adj   -- speed-corrected hGRF impulse
    
%% Set Constants

    g = 9.80665; %acceration due to gravity (m/s^2)
    
    F_v = m*g*(pi/2)*(t_f/t_c+1); %vGRF max estimate per Morin et al. (2005)
    d_y = (F_v*t_c^2)/(m*pi^2)-g/8*(t_c)^2; %vertical CoM disp. at midstance per Morin et al. (2005)
    
%% Numerically Solve for v_c

syms z

v_c = vpasolve(0 == (v-z*(t_c/(t_c+t_f)))*((t_c+t_f)/t_f)...
    -(1/(2*m*z))*(2*m*z^2+m*g*d_y-(m/8)*t_f^2*g^2+(F_v*(L-sqrt(L^2-(z*t_c/2)^2))^2+d_y)/(2*(L-sqrt(L^2-(z*t_c/2)^2)+d_y))), z, [0.8*v,v]);

v_c = double(v_c);

%% Solve for v_f

% Flight Velocity
v_f = (v-v_c*(t_c/(t_f+t_c)))*((t_c+t_f)/t_f);

%% Approximate F_h

% Peak hGRF (Braking and Propulsion)
F_h = (2*pi/t_c)*(v_f-v_c)*m;    % Absolute (N)
F_h_rel = F_h/(m*g);             % Relative (BW)

% hGRF Impulse
Imp_h = F_h*(t_c/pi);

%% Speed-Corrected hGRF

% Model Coefficients
v_center = 4.417989;      
vel_c = v_c - v_center;
B_0 = -0.04619409;
B_v = -0.06456883;

%Estimate-Actual
Fh_diff = B_0 + B_v*vel_c;

%Speed-corrected hGRF and Impulse
F_h_adj = F_h - Fh_diff*(m*g);
F_h_adj_rel = F_h_rel - Fh_diff;
Imp_h_adj = F_h_adj*(t_c/pi);

    end