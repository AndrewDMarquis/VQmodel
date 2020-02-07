function [] = COMPARE_NETWORK_PLOT(segL, COORD, Net, FSA, FSB, space, str, LFLAG, FN)


COORDA = COORD; COORDA(:,1) = COORD(:,1)-space;
COORDB = COORD; COORDB(:,1) = COORD(:,1)+space;

Nseg = size(segL,2); %number of segments
%%% make segment list for plotting
segP{Nseg} = []; %preallocate empty cell array to store "plot segment"
for i = 1:Nseg
   par    =  Net(segL{i},1);                     %find parents nodes
   TF     = logical(1 - ismember(par,segL{i}));  %logical index of parent node outside of segment
   c      = par(TF);                             %extract parent node index
   if c == 0
       segP{i} = segL{i}; %if c == 0, then this is the primary segment and we don't need to add the parent
   else
       segP{i} = [segL{i} c]; %"plot segment" is a cell array that has the parent node in the segment list - this ensure segments are connected together when plotted
   end
end

%%% shared color map
lowA = min(min(FSA));
hiA  = max(max(FSA));
lowB = min(min(FSB));
hiB  = max(max(FSB));
low = min([lowB,lowA]);
hi  = max([hiB,hiA]);

%%%color map
Nc = 100;      %number of colors
CM = jet(Nc);  %using the jet color map
m  = lowA;  %minimum flow
M  = hiA;  %maximum flow
dF = (M-m)/Nc; %partition size

CPA = ones(1,Nseg); %colors for each face;
CPB = ones(1,Nseg); %colors for each face;
for i = 1:Nc     %iterate by number of colors
    clear IND
    INDA = m+(i*dF) < FSA; %TF logical to dermine if color index needs to be increased
    INDB = m+(i*dF) < FSB; %TF logical to dermine if color index needs to be increased
    CPA(INDA) = CPA(INDA)+1; %increase color index if IND == 1
    CPB(INDB) = CPB(INDB)+1; %increase color index if IND == 1
end

if max([max(CPA) max(CPB)]) > Nc
   CM = jet(Nc+1); %sometimes we get an index +1 outside the OG colormap. If that happens, re-make the color map
end

%%% plot colorized network
figure(FN);
for ii = 1:Nseg
    hold on
    plot(COORDA(segP{ii},1),COORDA(segP{ii},2),'color', CM(CPA(ii),:),'linewidth',2)
    hold on
    plot(COORDB(segP{ii},1),COORDB(segP{ii},2),'color', CM(CPB(ii),:),'linewidth',2)
end

%%% colorbar
del = (hi - low)/10;
LAB{10} = [];
for x = 0:9
    if LFLAG == 1
        LAB{x+1} = sprintf('%0.2d',10^(low+x*del));
    else
        LAB{x+1} = round(low + x*del,2);
    end
end
colormap jet
c = colorbar('Location', 'southoutside');
set(c,'fontsize',18)
c.Label.String = str;
c.Ticks = (0:9)/9; %where to place tick marks
c.TickLabels = LAB;