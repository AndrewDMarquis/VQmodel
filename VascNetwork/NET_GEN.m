function [NetF, i3F, COORD, A] = NET_GEN(grow_size, K, Xc, Yc, root, display)
%%%% This function generates the topology and geometry and a vascular
%%%% network given a domain to cover, and some basic governing parameters
% Xc and Yc are vectors denoting the (x,y) coordiante pairs of seed points that cover the domain we want the network to span
% root is the coordinate pair denoting where we want to network to begin
% grow_size determines how far away a new node should be placed, and K is the number of new nodes the algorithm will generate
% display is an integer that controls how often to update the console with how far along the algorithm is - for practial use, display << K
%%%% last updated: 2/25/19 - Drew Marquis

% NetF is the primary output - it contains the spatial coordinates of all the nodes, the type of node, and the parent node each is connected to
% i3F is the list of indices corresponding to terminal nodes
% COORD is matrix containing the coordinates of all nodes - extracted from NetF
% A is an adjaceny matrix that represents the connectivey of the nodes - we
% treat the vascular tree as an undirected graph for ease of plotting with MATLAB's gplot function

%%%%% set up network data structure
%%%%% parent type   x                    y
Net = [0     1      root(1)              root(2)   
       1     3      root(1)-grow_size    root(2)]; %initialize network with first 2 nodes manually placed
%%% "parent" is the index of the node immediately connected (directed toward the root of the network) to the this node
%%% nodes can be type 1, 2, or 3. 1 - "normal node", 2 - bifurcation, 3 - terminal point of the network
%%% x and y are the spatial coordiantes for plotting

%%% network growing algorithm
Net = [Net; zeros(K,4)]; %stack root points and preallocate matrix to store node informatio
i3  = 2; % (initial) list of terminal nodes indices - updated through the following loop
for k = 1:K
    %%% find seed point with maximum distance to terminal nodes
    dist1           = sqrt((Xc - Net(i3,3)').^2 + (Yc - Net(i3,4)').^2); %distances between seed points and terminal nodes - organized as a matrix with (seed points) by (terminal node)
    [~,max_c_index] = max(min(dist1,[],2));                              %the operative search step - each terminal node has a seed point that it is furthest from. We want to the smallest of those furthest seed point distances
    
    %%% find node closest to (Xc,Yc)[max_c_index]
    dist2    = sqrt( (Xc(max_c_index) - Net(1:(k+1),3)).^2 + (Yc(max_c_index) - Net(1:(k+1),4)).^2 ); %distances between furthest point and ALL vascular network nodes
    [~,iNet] = min(dist2); %iNet is the index of closest node
    s        = [Xc(max_c_index) - Net(iNet,3) , Yc(max_c_index) - Net(iNet,4)];
    s        = s./norm(s); %unit vector pointing from terminal node to "new" terminal node.
    
    %%% we do different things based on what what type of node is closest 
    if Net(iNet,2) == 3  % nearest node is a terminal point
        Net(iNet,2) = 1; % reassign node type to normal 
        Net(k+2,:) = [iNet, 3, Net(iNet,3)+grow_size*s(1),  Net(iNet,4)+grow_size*s(2)]; %add new terminal node
    elseif Net(iNet,2) == 1 %nearest node is a normal node
        Net(iNet,2) = 2;    %reassign node type to bifurcation
        Net(k+2,:) = [iNet, 3, Net(iNet,3)+grow_size*s(1),  Net(iNet,4)+grow_size*s(2)]; %add new terminal
    elseif Net(iNet,2) == 2     %nearest node is a bifircation node
        while Net(iNet,2) ~= 1  %this algorithm is pretty stable - with the exception of this while loop. It is possible to go alllll the way back to the original node and '0' indices don't exist in matlab. Should only be an issue if number of seed points < number of nodes. Can get screwy based on the scaling of "grow_size" and domain of P 
            iNet = Net(iNet,1); %find a non-bifurcation parent node directed toward the root - we don't want to deal with bifurcations
        end 
        s = [Xc(max_c_index) - Net(iNet,3) , Yc(max_c_index) - Net(iNet,4)];
        s = s./norm(s);  %make new normal vector
        Net(iNet,2) = 2; %reassign node type to bifurcation
        Net(k+2,:) = [iNet, 3, Net(iNet,3)+grow_size*s(1),  Net(iNet,4)+grow_size*s(2)]; %add new terminal node
    end
    i3 = find(Net(:,2) == 3); %remake i3 - find all the terminal nodes
    
    %%% if k is a multiple of 'display', print the index in the console window. Nice to see if the
    %%% algorithm is making progress when it is taking a long time to run
    if mod(k,display) == 0
       disp(k) 
    end
end

%%% post-algorithm organization
%%% for ease of plotting, we treat the network as an undirected graph so we can use Matlab's gplot function
K      = K+2;        %total number of nodes
COORD  = Net(:,3:4); %coordinates stacked into a 2 by K matrix - extracted from Net
A      = zeros(K,K); %preallocated adjaceny matrix
for i = 2:K          %I bet it is possible to assemble this in the above loop, but conceptually I think this is cleaner
    A(i,Net(i)) = 1; %add connection to adjaceny matrix
end
A = A'+A; % yay symetric matrices. Depending on how many nodes you have and available RAM on your computer, it may be advantagoues to store A as a sparse matrix to conserve memory

%%% - "final" version of these variables to return - Net and i3 are updated throughout the algorithm and this ensures nothing is returned prematurely
NetF = Net;
i3F = i3;