%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = quarter_circle_domain(Demand,Arg)
BOXD = 9e3;
  BdBox = [-BOXD BOXD -BOXD BOXD];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)

%     s = BdBox(2);
    
    Aref = 6.714326703426993e+07;
    s = sqrt(4*Aref/pi);
    
    bc = dCircle(P,0,0,s);
    lr = dRectangle(P,0,s,-s,s);
    lr2 = dRectangle(P,-s,s,-s,0);
    
    Dist = dDiff(dDiff(bc,lr),lr2);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)

x = cell(2,1); % no boundary conditions
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%