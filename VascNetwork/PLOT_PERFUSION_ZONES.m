function [] = PLOT_PERFUSION_ZONES(Nt, MAP, Element, Node, Ft, str, LFLAG, fN)
%%%% this function plots the spatial perfusion zones associated with each
%%%% terminal branch. An important thing this function does is merge all
%%%% the seed points into one big polygon, so we don't have to plot the
%%%% individual polygon associated with each seed point
%%%  last updated: 2/25/19 - Drew Marquis

% new Element structure - rather than plotting a polygon for each seed
% point, we want to plot just the perfusion zone for each terminal branch
nElem{Nt} = []; 
for i = 1:Nt
    for j = 1:length(MAP{i})
        nElem{i} = [nElem{i} Element{MAP{i}(j)}]; %store indices of all nodes for a given terminal branch
    end
    [~,IND]  = unique(nElem{i});                        %find the unique indices
    nElem{i} = nElem{i}(sort(IND));                     %get rid of the repeats
    DT       = delaunayTriangulation(Node(nElem{i},:)); %graph theory is weird? unsure of what the DT struct does, but is needed for the next line
    k        = convexHull(DT);                          %the operative step - finds the exterior nodes
    nElem{i} = nElem{i}(k);                             %re-index with the exterior nodes - plot ready.
end

MaxNVer = max(cellfun(@numel,nElem));                   %max numumber of vertices in the entire mesh
PadWNaN = @(E) [E NaN(1,MaxNVer-numel(E))];             %Pad cells with NaN
ElemMat = cellfun(PadWNaN,nElem,'UniformOutput',false); %honesly not 100% sure how this and the previous line work - I copied/pasted this from PolyMesher source code
ElemMat = vertcat(ElemMat{:});                          %Create padded element matrix

%%%color map
Nc = 100;      %number of colors
CM = jet(Nc);  %using the jet color map
m  = min(Ft);  %minimum flow
M  = max(Ft);  %maximum flow
dF = (M-m)/Nc; %partition size

CP = ones(1,Nt); %colors for each face;
for i = 1:Nc     %iterate by number of colors
    clear IND
    IND = m+(i*dF) < Ft; %TF logical to dermine if color index needs to be increased
    CP(IND) = CP(IND)+1; %increase color index if IND == 1
end

if max(CP) > Nc
   CM = jet(Nc+1); %sometimes we get an index +1 outside the OG colormap. If that happens, re-make the color map
end

%%% actual plotting
figure(fN);
for k = 1:Nt
    clear temp
    temp = ElemMat(k,:); %extract the indices corresponding the verticies which define the perfusion zone
    hold on
    patch('Faces', temp, 'Vertices', Node, 'FaceColor', CM(CP(k),:), 'LineStyle', ':'); %plot the appropriately colored perfusion zone
end


low = min(min(Ft));
hi  = max(max(Ft));
if LFLAG == 1 %log scaled flag - use if you took the log10 of the Ft input
    low = 10^low; hi = 10^hi;
end
del = (hi - low)/10;
LAB{10} = [];
for x = 0:9
    if LFLAG == 0
        LAB{x+1} = round(low+x*del,1,'significant');
    else
        LAB{x+1} = sprintf('%0.2d',round(low + x*del,3,'significant'));
    end
end

colormap jet
c = colorbar('Location', 'eastoutside');
set(c,'fontsize',18)
c.Label.String = str;
c.Ticks = (0:9)/9; %where to place tick marks
c.TickLabels = LAB;
