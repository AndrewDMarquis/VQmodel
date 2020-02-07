function [] = COMPARE_3PERFUSION_ZONES(Nt, MAP, Element, Node, FtA, FtB, FtC, str, SPACE, LFLAG, colmap,fN)
%%% This function plots quantities related to the capillary compartments of
%%% VQ_RHS. Specifically for comparisons between 3 version of the same network

%make new node lists for each of the 3 things being compared
NodeA = Node; NodeA(:,1) = Node(:,1)-SPACE; 
NodeB = Node;
NodeC = Node; NodeC(:,1) = Node(:,1)+SPACE;

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

%%% shared color map
lowA = min(min(FtA));
hiA  = max(max(FtA));
lowB = min(min(FtB));
hiB  = max(max(FtB));
lowC = min(min(FtC));
hiC  = max(max(FtC));
low  = min([lowA,lowB,lowC]);
hi   = max([hiA,hiB,hiC]);

%%%colorbar bounds
if LFLAG == 1
    N = 10^abs(floor(low));
    low = log10(floor((10^low)*N)/N); %rounding down lowest value
    
    M = 10^(abs(floor(hi)));
    hi = log10(ceil(10^hi*M)/M); %rounding up highest value
else
    N = 10^abs(floor(log10(low)));
    low = floor(low*N)/N;
    
    M = 10^(abs(floor(log10(hi)))); %same as above but with non-log scaled inputs
    hi = ceil(hi*M)/M;
end

%%%color map
Nc = 1000;     %number of colors
m  = low;      %minimum flow
M  = hi;       %maximum flow
dF = (M-m)/Nc; %partition size

if colmap == 1 %choose what predefined colormap you want to use. If you change it here, also change the colorbar option further down
    CM = hot(Nc);
elseif colmap == 2
    CM = pink(Nc);
elseif colmap == 3
    CM = jet(Nc);
end

CPA = ones(1,Nt); %colors for each face;
CPB = ones(1,Nt); %colors for each face;
CPC = ones(1,Nt); %colors for each face;
for i = 1:Nc      %iterate by number of colors
    clear IND
    INDA = m+(i*dF) < FtA;   %TF logical to dermine if color index needs to be increased
    INDB = m+(i*dF) < FtB;   %TF logical to dermine if color index needs to be increased
    INDC = m+(i*dF) < FtC;   %TF logical to dermine if color index needs to be increased
    CPA(INDA) = CPA(INDA)+1; %increase color index if IND == 1
    CPB(INDB) = CPB(INDB)+1; %increase color index if IND == 1
    CPC(INDC) = CPC(INDC)+1; %increase color index if IND == 1
end

if max([max(CPA) max(CPB) max(CPC)]) > Nc
   CM = jet(Nc+1); %sometimes we get an index +1 outside the OG colormap. If that happens, re-make the color map a little bigger
end

%%% actual plotting
figure(fN);
for k = 1:Nt %iterate by perfusion zone
    clear temp
    temp = ElemMat(k,:); %extract the indices corresponding the verticies which define the perfusion zone
    hold on
    patch('Faces', temp, 'Vertices', NodeA, 'FaceColor', CM(CPA(k),:), 'LineStyle', ':'); %plot the appropriately colored perfusion zone - input A
    patch('Faces', temp, 'Vertices', NodeB, 'FaceColor', CM(CPB(k),:), 'LineStyle', ':'); %plot the appropriately colored perfusion zone - input B
    patch('Faces', temp, 'Vertices', NodeC, 'FaceColor', CM(CPC(k),:), 'LineStyle', ':'); %plot the appropriately colored perfusion zone - input C
end


%%% color bar labels
if LFLAG == 1
    [ticks,tickLabels,scaleStr,minorTicks] = calcticks(exp([low hi]),'y',1,'log'); %calcticks is neat function I found on mathworks file exchange to make tick labels
else
    [ticks,tickLabels,scaleStr,minorTicks] = calcticks([low hi]);
end

if colmap == 1
    colormap hot
elseif colmap == 2
    colormap pink
elseif colmap == 3
    colormap jet
end

c = colorbar('Location', 'westoutside');
set(c,'fontsize',22)
c.Label.String = str;

if length(tickLabels) < 3
    c.TickLabels = num2str(minorTicks',2);
else
    c.TickLabels = tickLabels; %applying new tick labels
end
