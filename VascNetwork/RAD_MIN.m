function OUT = RAD_MIN(r,L,F,dP)
%%% this is a cost function used to determine the radius of branches in the
%%% fractal tree. Since we're going to eventually use an non-linear relationship
%%% to define viscosity as a function of vessel diamter, we want to account
%%% for that now (see VISCOSITY function). All of this is to say we're
%%% going use a root finder to find what radius makes OUT = 0.

H = 0.45; % 45% hematocrit

eta_p = 1.2;               %viscosity of plasma (mPa*s) - measured value from: Késmárky, Gábor, et al. "Plasma viscosity: a forgotten variable." Clinical hemorheology and microcirculation 39.1–4 (2008): 243-246.
eta_r = VISCOSITY(2.*r,H); %relative viscosity (unitless) as a function of diameter and hematocrit
eta   = eta_r*eta_p*1e-3;  %apparent viscostiy of blood (Pa*s)

dP    = dP*133.322;        %converting mmHg to Pa (1 mmHg = 133.322 Pa)

temp  = (8*eta*L*F)./(pi*dP); %algebra to solve for r from Poisuelle's law
OUT   = nthroot(temp,4)-r;    %fzero to find the value of r that makes OUT = 0
