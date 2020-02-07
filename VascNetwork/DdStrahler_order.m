function [SOF, DOF, iterF] = DdStrahler_order(Net, segL, R)
%%% This function determines the diamter definied Strahler ordering of the segments in the
%%% vascular network made by NET_GEN. Note that this function requires you to determine the 
%%% radius/diamter of each segment via the ANATOMY_STATS function before calling this one
%%% last updated: 2/25/19 - Drew Marquis

%Net  - primary data strucutre from NET_GEN
%segL - cell array were each entry is a segment with a vector of inidcies corresponding to the nodes in Net - made in ANATOMY_STATS
%R    - radius of each segment - made in ANATOMY_STATS

%SOF   - Strahler Order of each segment
%DOF   - diameters in each Stahler order
%iterF - number of iterations we need to use the diamter definied rules to reach a converged Stahler Ordering for the entire network

Ns = length(R);                   %number of segments
Nt = length(find(Net(:,2) == 3)); %number of terminal segments
D  = 2*R;                         %diameters of segments

%%% Strahler Ordering
SO       = zeros(1,Ns); %preallocate vector to store strahler orders for each segment- indices correspond to the cell-array entires of segL and indices of R
SO(1:Nt) = 1;           %terminal segments are all order 1
count    = Nt;          %we've assigned an order for Nt vessels, so our count is Nt
IND      = (1:Nt)';     %vector of indices corresponding to current vessel segments
while count < Ns
    clear i k m BIF
    %%% find bifurcation points
    BIF = zeros(1,length(IND)); %preallocate list of bifucation indices
    for i = 1:length(IND)       %iterate by terminal nodes
        clear par TF c j
        par    = Net(segL{IND(i)},1);                      %list of parent nodes
        TF     = logical(1 - ismember(par,segL{IND(i)}));  %logical index of parent node outside of segment
        c      = par(TF);                                  %extract parent node index
        for j = 1:Ns
            if ismember(c,segL{j}) %loop to determine which segment that parent is in
                BIF(i) = j;        %store bifurcation index
                break;
            end
        end
    end
    %%% find the bifurcations we can actually handle this iteration
    ubif = []; %preallocate empty array to store useful BIF indices
    dau  = []; %empty array to store daughter segment indices
    for k = 1:length(BIF)
        TF2 = BIF == BIF(k); %compare entires of BIF with itself
        if sum(TF2) == 2
            ubif = [ubif BIF(k)];  %store useful bif index
            dau  = [dau IND(TF2)]; %store coresponding daughter index
        end 
    end
    %%% determine Strahler order (NOT diameter defined)
    [ubif,ind] = unique(ubif);            %by design, we double count the useful bifurcation nodes - so we need to use 'unique' to remove redundancy from the list(s)
    dau        = dau(:,ind);              %remove redundant entries by indexing with output from previous line
    for m = 1:length(ubif)                %iterate by useful bifurcation points
        if SO(dau(1,m)) == SO(dau(2,m))   %if orders of daughters are equal
            SO(ubif(m)) = SO(dau(1,m))+1; %increase order of confluent parent by 1
        else
            SO(ubif(m)) = max(SO(dau(:,m))); %if orders of daughters are not the same, confluent parent order is the maximum of the daughter vessels
        end
        ti1 = IND == dau(1,m); %index of IND to update with new parent
        ti2 = IND == dau(2,m); %index of IND to delete - the choice between ti1 and and ti2 is arbitrary and can be used interchangeably
        IND(ti1) = ubif(m);    %updating a node -  which node is updated and which is deleted is arbitrary, but it is important to update a node BEFORE you delete the other, otherwise indexing may be offset and that is no good
        IND(ti2) = [];         %deleting the other node - see the comment in the line above
        count    = count+1;    %update counter
    end
end

%%% diamter definied (dd) Strahler Ordering
SOdd = zeros(1,Ns); %preallocate dd Stahler orders

MO     = max(SO);     %find biggest order - "Max Order"
DO{MO} = [];          %cell array to store diameters for each order
Dmean  = zeros(1,MO); %preallocate mean diameter for each order
Dsd    = zeros(1,MO); %preallocate standard deviation for each order

change = 1; %change in ordering as we use dd rules
iter   = 0; %number of times we have to iterate with dd rules
while change > 0.01 %we say that the ordering scheme has converged if there is less than 1% change in the orders
    clear i
    for i = 1:MO                %iterate by strahler orders
        clear ind
        ind      = SO == i;     %indices corressponding to i'th order
        DO{i}    = D(ind);      %find the diameters and store in cell array
        Dmean(i) = mean(DO{i}); %mean
        Dsd(i)   = std(DO{i});  %standard deviation
    end
    DB = (Dmean(2:end)+Dsd(2:end) + Dmean(1:end-1) - Dsd(1:end-1))/2;  %diameter boundaries to partition orders
    nB = Dmean(end) + Dsd(end); %'new Boundary' - in the event that re-ordering produces orders bigger than MO
    clear j
    for j = 1:Ns                %iterate by each segment
        clear TF nO
        if D(j) > nB %if j'th segment diamter is bigger than the largest boundary of the current largest Strahler order
            MO      = MO+1;            %increase max Strahler order
            SOdd(j) = MO;              %assign dd Strahler order
            DB      = [DB nB];         %stack new boundary with diamter boundaries
            nB      = D(j) + Dsd(end); %make a bigger nB to account for another possible new order - this is super crude and can probably be improved? 10/18/18
        else
            TF      = D(j) < DB;    %compare j'th diameter to bounds
            nO      = MO - sum(TF); %mapping to new order
            SOdd(j) = nO;           %assign dd Strahler order
        end
    end
    x      = find(SOdd - SO); %indices of segments who's Strahler order changed
    change = length(x)/Ns;    %percent change (number of changed segments/total number of segments)
    SO     = SOdd;            %assign SO as SOdd so we can re-iterate through the while loop
    iter   = iter + 1;        %add 1 to iter.
    disp(['ddSO: ', num2str(iter), ' ', num2str(change)])
end

%%%% re-making DO
clear i
for i = 1:MO                %iterate by strahler orders
    clear ind
    ind      = SOdd == i;   %indices corressponding to i'th order
    DO{i}    = D(ind);      %find the diameters and store in cell array
end

%%% sometimes we get empty 1rst and 2nd orders (potentially more) after the diamter defined upward shuffling.
%%% This gets rid of those empty cell arrays
clear i ind
ind = [];
for i = 1:MO %find the largest empty strahler order
    if isempty(DO{i})
        ind = i; 
    end
end
if ~isempty(ind) %only execute this if we have empty orders - note that this logical condition is determining if ind is the empty array we preallocated
    nDO{MO-ind} = [];
    for j = MO:-1:(ind+1)
        nDO{j-ind} = DO{j}; %recollecting the orders
    end
    DO   = nDO;
    SOdd = SOdd-ind;
end

%%%% 'final' version to prevent premature return of outputs
SOF   = SOdd; 
DOF   = DO;
iterF = iter;