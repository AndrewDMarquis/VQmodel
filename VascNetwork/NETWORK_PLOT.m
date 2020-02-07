function [] = NETWORK_PLOT(segL, COORD, Net, FS, str, LFLAG, FN)
%%% This function plots a colorized version of the vascular network. Useful
%%% for visualizing flow/pressure/variable distribution in large
%%% networks. The colormap is based on the values in the vector FS
%%% corresponding to the segmetns in segL. Entries in segL
%%% correspond the node coordinates in the (x,y) plane to plot the network.
%%% Net is the data struct cotaining infomration on how network nodes are
%%% connected to one another
%%% last updated: 2/25/19 - Drew Marquis

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
%%% make colormap and assign colors to each segment
Nc = 100;      %number of colors
CM = jet(Nc); %using color map
m  = min(FS);  %minimum flow
M  = max(FS);  %maximum flow
dF = (M-m)/Nc; %partition size

CP = ones(1,Nseg); %colors for each vessel segment
for j = 1:Nc     %iterate by number of colors
    clear IND
    IND = m+(j*dF) < FS; %TF logical to dermine if color index needs to be increased
    CP(IND) = CP(IND)+1; %increase color index if IND == 1
end
if max(CP) > Nc
   CM = jet(Nc+1);
end

hi = max(FS);
low = min(FS);
%%% plot colorized network
figure(FN);
for ii = 1:Nseg
    hold on
    plot(COORD(segP{ii},1),COORD(segP{ii},2),'color', CM(CP(ii),:),'linewidth',2)
end
axis equal

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