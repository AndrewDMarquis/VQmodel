function eta = VISCOSITY(D,Hd)
%%% This function defines the rheological properies of blood viscosty accounting for
%%% the Fahraeus and Fahraeus-Lendquvist effects
%%% Output the relative viscostiy of blood (unitless) as function of vessel diamter
%%% (um) and hematocrit (%)
%%% This empirical relationship was derived and fit to experimetnal data in
%%% "Resistance to Blood Flow in Microvessels In Vivo" - A.R.Preis et al.(1994) AHA journal
%%% also referred to as the "in vivo viscosity relationship"

C        = (0.8+exp(-075*D)).*(-1+1./(1+1e-11*D.^12))+1./(1+1e-11*D.^12);
eta_star = 6*exp(-0.085.*D)+3.2-2.44.*exp(-0.06*D.^0.645);
eta      = (1+(eta_star-1).*(((1-Hd).^C-1)./(((1-0.45).^C)-1)).*(D./(D-1.1)).^2).*(D./(D-1.1)).^2;


