est.vc <- function(t_c, t_f, v, m, L){

# Author:
# Geoffrey Burns
# University of Michigan
  
# Inputs:
# t_f -- flight time (s)
# t_c -- contact time (s)
# v   -- subject's average running speed (m/s)
# m   -- subject's mass (kg)
# L   -- subject's leg length (m)
#        Approximate leg length via height: L = 0.53*h
#        where h is the subject's height (m) per Winter (2005)

# Outputs:
# F_h     -- peak hGRF (N) (negative in braking/positive in propulsion)
# F_h.rel -- peak hGRF (BW) (negative in braking/positive in propulsion)
# Imp_h   -- hGRF impulse (N-s) (negative in braking/positive in propulsion)
# F_h.adj     -- corrected peak hGRF (N) (negative in braking/positive in propulsion)
# F_h.rel.adj -- corrected peak hGRF (BW) (negative in braking/positive in propulsion)
# Imp_h.adj   -- corrected hGRF impulse (N-s) (negative in braking/positive in propulsion)

    
## Set Constants

    g <- 9.80665; #acceration due to gravity (m/s^2)
    
    F_v <- m*g*(pi/2)*(t_f/t_c+1); #vGRF max estimate per Morin et al. (2005)
    d_y <- (F_v*t_c^2)/(m*pi^2)-g/8*(t_c)^2; #vertical CoM disp. at midstance per Morin et al. (2005)
    
## Numerically Solve for v_c

f <- function(z){
  (v-z*(t_c/(t_c+t_f)))*((t_c+t_f)/t_f)-(1/(2*m*z))*(2*m*z^2+m*g*d_y-(m/8)*t_f^2*g^2+(F_v*(L-sqrt(L^2-(z*t_c/2)^2))^2+d_y)/(2*(L-sqrt(L^2-(z*t_c/2)^2)+d_y)))
}

fvec <- Vectorize(f)

v_c.solve <- uniroot(f, c(0.8*v,v));

v_c <- v_c.solve$root

## Solve for v_f and v_ms

# Flight Velocity    
v_f <- (v-v_c*(t_c/(t_f+t_c)))*((t_c+t_f)/t_f);

## Approximate F_h

# Peak hGRF (Braking and Propulsion)
F_h <- (2*pi/t_c)*(v_f-v_c)*m # Absolute (N)
F_h.rel <- F_h/(m*g)    # Relative (BW)
# hGRF Impulse
Imp_h <- F_h*(t_c/pi);

## Adjusted hGRF

#Model Coefficients
v.center <- 4.417989
vel.c <- v_c - v.center

B.0 <- -0.04619409
B.v <- -0.06456883

#Estimate-True
Fh.diff <- B.0 + B.v*vel.c

#Speed-corrected hGRF and Impulse
F_h.adj <- F_h - Fh.diff*(m*g)
F_h.adj.rel <- F_h.rel - Fh.diff
Imp_h.adj <- F_h.adj*(t_c/pi)

## Results List

output <- list(F_h = F_h,
               F_h.rel = F_h.rel,
               Imp_h = Imp_h,
               F_h.adj = F_h.adj,
               F_h.adj.rel = F_h.adj.rel,
               Imp_h.adj = Imp_h.adj)
output
}