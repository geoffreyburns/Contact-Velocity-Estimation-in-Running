function [v_c, v_f, v_ms] = vc_est(t_c, t_f, v, m, h)

%% Stance Velocity Estimation

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
% v_c       -- average horz. speed during contact (m/s)
% v_f       -- average horz. speed during flight (m/s)
% v_ms      -- minimum horz. speed at midstance (m/s)
    
%% Set Constants
    
    g = 9.80665; %acceration due to gravity (m/s^2)
    
    F_v = m*g*(pi/2)*(t_f/t_c+1); %vGRF max estimate per Morin et al. (2005)
    d_y = (F_v*t_c^2)/(m*pi^2)-g/8*(t_c)^2; %vertical CoM disp. at midstance per Morin et al. (2005)
    
%% Numerically Solve for v_c

syms z

v_c = vpasolve(0 == (v-z*(t_c/(t_c+t_f)))*((t_c+t_f)/t_f)...
    -(1/(2*m*z))*(2*m*z^2+m*g*d_y-(m/8)*t_f^2*g^2+(F_v*(L-sqrt(L^2-(z*t_c/2)^2))^2+d_y)/(2*(L-sqrt(L^2-(z*t_c/2)^2)+d_y))), z, [0.8*v,v]);

v_c = double(v_c);

%% Solve for v_f and v_ms

% Flight Velocity
v_f = (v-v_c*(t_c/(t_f+t_c)))*((t_c+t_f)/t_f);
% Midstance Velocity
v_ms = 2*v_c-v_f;

    end