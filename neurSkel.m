classdef neurSkel
    properties
        nodes % struct array, each node has the field 'links' (array of indices of all
              % other nodes connected to this node) and 'coords' (x, y, z
              % coordinates of this noe)
        links % n x 2 matrix, n = number of nodes. columns are the node indices at the end of
              % each link. Note some links are [0 0] due to some nodes having only one link
        linkDists % n x 2 matrix, column 1 is real distance, column 2 is electrotonic distance
        skel % stores the original table from the skeleton .csv file provided to the constructor
        synSet % a struct array
        % fields of synSet (stands for "synapse set"):
            % uniquePartners: holds the body_Ids of all unique partners in
            %   this synSet
            % uniquePartnerIndices
            % synLocs 
            % synLinks
            % synGroupings
            % rawData
            % synROIs
            % synROIkey
            % note: synLocs, synLinks apply to UNIQUE synapse locations.
            % Each location probably has many synapses (given by
            % synGroupings)
        
        synSetNames % human-understandable annotations for what each synSet is
        filepath % name of original .csv file given to the constructor
        synGraphs % cell array, each cell corresponds to one synSet
                  % each synSet's synGraph is an nx1 cell array where n is
                  % the number of unique synapse locations in the synSet.
                  % each cell of that nx1 cell array is an mx3 matrix where
                  % m is the number of "synapse neighbors" that synapse
                  % location has.
                  % so obj.synGraphs{1}{302} is synapse location 302 in
                  % synSet 1, and might have, say, 4 synapse neighbors if
                  % it's close to a branching point. 
                  % the 3 columns are: 
                  % 1. the index of the neighboring synapse location
                  % 2. the actual distance along the skeleton
                  % 3. the 'electrotonic' distance (distance/sqrt(radius) for cylinders;
                  %    2*distance*(sqrt(radius1)-sqrt(radius2))/(radius1-radius2) for truncated cones)
        labeledNodes % label nodes according to if they belong to disconnected fragments
                     % e.g. if a skeleton has 1 big continuous branch but 3
                     % fragments, each fragment would have a number 1, 2, 3
                     % or 4
        labeledLinks % same as for labeledNodes, but for links
        directions % 1 = this neuron is pre; 2 = this neuron is post. This is an nx1 array like synSetNames, 
                   % the ith number is the direction of the ith synSet
                   
        nodesInROI
        linksInROI
    end
    
    methods
        function obj = neurSkel(filepath)
            % constructor method takes in the file name of a .csv file
            % storing the entire skeleton of a neuron
            if nargin>0
                obj.filepath = filepath;
                obj.skel = readtable(obj.filepath);
                obj = obj.getNodes();
            end
        end
        
        function obj = getNodes(obj)
            % create nodes and links from the raw skeleton table
            disp('creating nodes...');
            tic
            numNodes = size(obj.skel,1);
            % preallocate nodes
            obj.nodes = struct('coords',cell(numNodes,1),'links',cell(numNodes,1));
            obj.links = zeros(numNodes,2);
            tic
            for i=1:numNodes
                % store the coordinates of this node
                obj.nodes(i).coords=[obj.skel.x(i) obj.skel.y(i) obj.skel.z(i)];
            end
            disp('creating links...');
            for i=1:numNodes
                if mod(i,10000)==0
                    i
                end
                %store the link between this node and its neighbour (-1 means no link)
                if (obj.skel.link(i)>0)
                    obj.nodes(i).links=[obj.nodes(i).links obj.skel.link(i)];
                    obj.nodes(obj.skel.link(i)).links = [obj.nodes(obj.skel.link(i)).links i];
                    % add the link to links
                    obj.links(i,:) = [i obj.skel.link(i)];
                    % calculate distances
                    % real distance
                    obj.linkDists(i,1) = pdist2(obj.nodes(obj.links(i,1)).coords, ...
                            obj.nodes(obj.links(i,2)).coords);
                    % electrotonic distance
                    radius1 = obj.skel.radius(obj.links(i,1));
                    radius2 = obj.skel.radius(obj.links(i,2));
                    if radius1 ~= radius2
                        obj.linkDists(i,2) = obj.linkDists(i,1)*2*(sqrt(radius1)-sqrt(radius2))/(radius1-radius2);
                    else
                        obj.linkDists(i,2) = obj.linkDists(i,1)/sqrt(radius1);
                    end
                end
            end
            toc
        end
        
        function obj = labelNodes(obj)
            % tag links according to if they are disconnected from the rest
            % of the skeleton
            numNodes = length(obj.nodes);
            Q = 1:numNodes;
            obj.labeledNodes = zeros(numNodes,1);
            obj.labeledLinks = zeros(size(obj.links,1),1);
            index = 1;
            visited = zeros(numNodes,1);
            % use a queue to walk through all connected nodes
            while ~isempty(Q)
                % if mod(length(Q),100)==0
                %     fprintf("Q now has %d elements", length(Q));
                % end
                % disp('length(Q)=%d\n',length(Q));
                stack = Q(1);
                while ~isempty(stack)
                    node = stack(1);
                    stack = stack(2:end);
                    if ~visited(node)
                        visited(node) = 1;
                        obj.labeledNodes(node) = index;
                        obj.labeledLinks(obj.links(:,1)==node) = index;
                        Q = Q(Q~=node);
                        toPush = obj.nodes(node).links;
                        stack = [obj.nodes(node).links, stack];
                    end
                    if mod(nnz(visited),1000)==0
                        fprintf("visited %d nodes\n", nnz(visited));
                    end
                end
                index = index+1;
            end
        end
        
        function visited = traverseSynGraph(obj, i, start)
            % for debugging: test if for the i-th synGraph, if you start
            % from synapse with the index given by "start", can you reach
            % every other synapse location in the i-th synSet by walking
            % along the neuron's skeleton? If not, either the synapses are
            % on disconnected fragments, or something is wrong with the
            % function getSynSetGraph. Note: the constant
            % "distanceTolerance" needs to be set correctly to prevent
            % errors
            graph = obj.synGraphs{i};
            numNodes = length(graph);
            stack = start;
            visited = false(numNodes,1);
            f = waitbar(0,'traversing...');
            while ~isempty(stack)
                node = stack(1);
                if mod(nnz(visited),1000)==0
                    waitbar(nnz(visited)/numNodes,f);
                end
                stack = stack(2:end);
                if ~visited(node)
                    visited(node) = true;
                    if ~isempty(graph{node})
                        toPush = graph{node}(:,1);
                        stack = [toPush; stack];
                    end
                end
            end
            close(f);
        end
        
        function result = findAsymSynNeighbors(obj, i)
            % this function is also for debugging getSynSetGraph
            graph = obj.synGraphs{i};
            result = zeros(length(graph),1);
            for j=1:length(graph)
                % for all the neighbors of j
                for k=1:size(graph{j},1)
                    % does k also have j as a neighbor?
                    kthneighbor = graph{j}(k,1);
                    if ~sum(graph{kthneighbor}(:,1)==j)
                        result(j) = 1;
                    end
                end
            end
        end

        function obj = removeStrayFragments(obj)
            % remove all disconnected fragments other than the longest one
            % not used in practice
            if isempty(obj.labeledNodes)||isempty(obj.labeledLinks)
                error('run neurSkel.labelNodes first!');
            end
            if size(obj.labeledNodes)~=size(obj.labeledLinks)
                error('labeledNodes and labeledLinks sizes do not match');
            end
            if mode(obj.labeledNodes)~=mode(obj.labeledLinks)
                error('modes do not match');
            end
            if isempty(obj.synSet)
                error('you have to localize the synapses to links before removing stray fragments');
            end
            % back up your object in case this removal goes horribly
            % wrong!!
            save(strcat('neurSkelbackup',datestr(now,'yyyymmddTHHMMSS'),'.mat'),'obj');
            modeLabel = mode(obj.labeledNodes);
            nodesToRemove = (obj.labeledNodes~=modeLabel);
            linksToRemove = (obj.labeledLinks~=modeLabel);
            obj.nodes(nodesToRemove) = [];
            obj.links(linksToRemove,:) = [];
            % retain labeledNodes and labeledLinks as 'relics' to remind us
            % of what the original indices were in the raw data
            obj.skel(nodesToRemove,:) = [];
            for i=1:length(obj.synSet)
                [~,synsToRemove,~] = intersect(obj.synSet(i).synLinks, find(linksToRemove));
                obj.synSet(i).synGroupings(synsToRemove) = [];
                obj.synSet(i).synLocs(synsToRemove,:) = [];
                obj.synSet(i).synLinks(synsToRemove) = [];
                % don't touch uniquePartnerIndices because synGroupings
                % still refers to the original indices in uniquePartnerIndices
            end
        end
        
        function obj = getSynSets(obj, filepaths, options)
            % This function creates the field 'synSet'
            % Inputs:
            % filepath: cell array of .csv file names
            % options: array defining the direction of connectivity for
            % each .csv file
            % option is 1 for the skeleton neuron = pre
            % option is 2 for the skeleton neuron = post
            % It maps every synapse to a location on the neuron's skeleton,
            % then groups together synapses with the same location into a
            % "synLoc" (i.e. unique synapse location). NB this includes
            % both genuinely polyadic synapses and also separate synapses
            % that just map to the same location
            % The function calls "findSynLocs"
            [rootpath,rootname,ext] = fileparts(obj.filepath);
            obj.synSet = struct;
            for i=1:length(filepaths)
                % say root is 'APL.csv'
                % and this set of synapses is 'APLtoKC'
                % then we want the files to be 'APLisPre-APLtoKC-synLocs.csv', 'APLisPre-APLtoKC-synLinks.csv',
                [~,name,~] = fileparts(filepaths{i});
                synLocsName = strcat(name,'-synLocs');
                synLinksName = strcat(name,'-synLinks');
                switch options(i)
                    case 1
                        directionTag = 'isPre';
                    case 2
                        directionTag = 'isPost';
                    otherwise
                        error('invalid options tag');
                end
                synLocsPath = strcat(rootpath,rootname,directionTag,'-',synLocsName,ext);
                synLinksPath = strcat(rootpath,rootname,directionTag,'-',synLinksName,ext);
                if 0 %isfile(synLocsPath)&&isfile(synLinksPath)
                    % left-over code for reading in pre-calculated
                    % locations, because it takes so long to run, but in
                    % practice not necessary, just save the neurSkel object
                    % in a .mat file
                    obj.synSet(i).synLocs = csvread(synLocsPath);
                    obj.synSet(i).synLinks = csvread(synLinksPath);
                else
                    disp(['getting synapse locations from ',filepaths{i}]);
                    synapses = readtable(filepaths{i});
                    % This is a hack for backward compatibility
                    if ismember('x_pre', synapses.Properties.VariableNames)
                        synapses.Properties.VariableNames{'x_pre'} = 'ux';
                        synapses.Properties.VariableNames{'y_pre'} = 'uy';
                        synapses.Properties.VariableNames{'z_pre'} = 'uz';
                        synapses.Properties.VariableNames{'x_post'} = 'dx';
                        synapses.Properties.VariableNames{'y_post'} = 'dy';
                        synapses.Properties.VariableNames{'z_post'} = 'dz';
                        synapses.Properties.VariableNames{'bodyId_pre'} = 'upstream_bodyId';
                        synapses.Properties.VariableNames{'bodyId_post'} = 'downstream_bodyId';
                        synapses.Properties.VariableNames{'roi_pre'} = 'upstream_roi';
                        synapses.Properties.VariableNames{'roi_post'} = 'downstream_roi';
                    end
                    switch options(i)
                        case 1 % if skeleton neuron is pre, we want the upstream partner
                            synOrigLocs = synapses{:,{'ux','uy','uz'}};%[synapses.ux(i) synapses.uy(i) synapses.uz(i)];
                            [obj.synSet(i).uniquePartners, ~, obj.synSet(i).uniquePartnerIndices] ...
                                = unique(synapses{:,{'downstream_bodyId'}});
                            % uniquePartners(uniquePartnerIndices) =
                            % original raw data
                        case 2 % if skeleton neuron is post, we want the downstream partner
                            synOrigLocs = synapses{:,{'dx','dy','dz'}};%[synapses.dx(i) synapses.dy(i) synapses.dz(i)];
                            [obj.synSet(i).uniquePartners, ~, obj.synSet(i).uniquePartnerIndices] ...
                                = unique(synapses{:,{'upstream_bodyId'}});
                        otherwise
                            error('parameter "option" needs to be 1 or 2');
                    end
                    obj.synSet(i).rawData = synapses;
                    
                    % group synapses that are in the same place (polyadic synapses)
                    [uniqueLocs, ~, ic] = unique(synOrigLocs, 'rows', 'stable');
                    % according to documentation of 'unique':
                    % synOrigLocs = uniqueLocs(IC,:)
                    obj.synSet(i).synGroupings = cell(length(uniqueLocs),1);
                    % each cell in synGroupings will get a list of indices in the
                    % original raw data that have a synapse in this same location
                    for j=1:length(uniqueLocs)
                        obj.synSet(i).synGroupings{j} = find(ic==j);
                    end
                    
                    %[obj.synSet(i).synLocs, obj.synSet(i).synLinks] = obj.findSynLocs(uniqueLocs);
                    [indivLocs, indivLinks] = obj.findSynLocs(synOrigLocs);
                    % group synapses that are in the same place 
                    % NOTE: this includes BOTH genuinely polyadic synapses
                    % (a single pre-synaptic density targeting multiple
                    % post-synaptic neurons) AND two separate synapses that
                    % happen to map to the same location on the skeleton
                    [obj.synSet(i).synLocs, ia, ic] = unique(indivLocs, 'rows', 'stable');
                    % according to documentation of 'unique':
                    % [C,IA,IC] = unique(A) also returns index vectors IA and IC such that
                    %     C = A(IA) and A = C(IC)
                    % indivLocs = uniqueSynLocs(ic,:)
                    % uniqueSynLocs = indivLocs(ia,:)
                    obj.synSet(i).synLinks = indivLinks(ia);
                    
                    obj.synSet(i).synGroupings = cell(length(obj.synSet(i).synLinks),1);
                    % each cell in synGroupings will get a list of indices in the
                    % original raw data that have a synapse in this same location
                    for j=1:length(obj.synSet(i).synLinks)
                        obj.synSet(i).synGroupings{j} = find(ic==j);
                    end
                    
                    %csvwrite(synLocsPath, obj.synSet(i).synLocs);
                    %csvwrite(synLinksPath, obj.synSet(i).synLinks);
                end
                obj.synSetNames{i} = strcat(rootname,directionTag,'-',name);
                obj.directions = options;
            end
        end
        
        function obj = labelSynsByROI(obj, synSetIndices)
            % labels every synapse with the ROI that it's tagged with in
            % the original .csv file 
            % this should be unnecessary - could have read the ROIs of the
            % partner neuron to start with in getSynSet! 
            is = synSetIndices;
            for i=1:length(is)
                switch obj.directions(is(i))
                    % there are occasional synapses that lie right at the
                    % boundary of 2 ROIs: pre is one ROI and post is in another
                    case 1
                        [obj.synSet(is(i)).synROIskey, ~, origSynROIs] = ...
                            unique(obj.synSet(is(i)).rawData{:,'downstream_roi'});
                    case 2
                        [obj.synSet(is(i)).synROIskey, ~, origSynROIs] = ...
                            unique(obj.synSet(is(i)).rawData{:,'upstream_roi'});
                end
                obj.synSet(is(i)).synROIs = zeros(length(obj.synSet(is(i)).synLinks),1);
                for j=1:length(obj.synSet(is(i)).synLinks)
                    roiLabels = origSynROIs(obj.synSet(is(i)).synGroupings{j});
                    % take the mode because the different synapses might
                    % not all be in the same ROI:
                    % this could happen if a neuron synapses onto
                    % multiple neurons from the same presynaptic
                    % density, but the synapse is on an ROI boundary
                    % and one target is in one ROI, other is in another
                    obj.synSet(is(i)).synROIs(j) = mode(roiLabels);
                end
            end
        end
        
        
        function [synLocs, synLinks] = findSynLocs(obj, synOrigLocs)
        % take in a list of 3d coordinates and find all their locations on the
        % skeleton
        % returns the 3d coordinates of where they map to on the skeleton,
        % and the index of the link on the skeleton on which the synLocs lie 

            numSyn = size(synOrigLocs,1);
            synLocs = zeros(numSyn,3);
            synLinks = zeros(numSyn,1);
            
            % ignore links that have 0 values
            linkIndices = 1:size(obj.links,1);
            nonzeroLinkIndices = linkIndices(all(obj.links,2));
            actualLinks = obj.links(nonzeroLinkIndices,:);
            
            % this is a vectorized version of the code commented out in the for loop.
            % this part doesn't need to be in the loop because A,B,N are fixed for the
            % skeleton, don't vary across which synapse is being tagged.
            % A is the points on the "first" end of the link
            A = obj.skel{actualLinks(:,1),{'x','y','z'}};
            % B is the points on the other end of the link
            B = obj.skel{actualLinks(:,2),{'x','y','z'}};
            % N is the vector pointing from A to the other end of the
            % link
            N = B - A;
            normN = sqrt(sum(N.^2,2)); % take the norm of each row of N
            unitN = N./repmat(normN,1,size(N,2));
            % tic
            for i=1:numSyn
                if mod(i,1000)==0
                    i
                end
                
                P = synOrigLocs(i,:);
                
                % dot product along dimension 2 (the vectors being dotted together are rows)
                dotProducts = min(normN,max(0,dot(P-A,unitN,2)));
                closestPointPerLink = A + repmat(dotProducts,1,3).*unitN;
                % note: this index refers to the index within actualLinks, not the
                % original links!
                [~,closestNonzeroLink] = min(pdist2(closestPointPerLink,P));
                % this line re-indexs back into the original link indices
                closestLink = nonzeroLinkIndices(closestNonzeroLink);
                synLocs(i,:) = closestPointPerLink(closestNonzeroLink,:);
                synLinks(i) = closestLink;
                % this code is vectorized above - the for loop runs 100x slower!!
                %     for k = nonzeroLinkIndices
                %
                %         % the point on one end of the link
                %         A = nodes(links(k,1)).coords;
                %         % heuristic to speed up processing: only search for links whose
                %         % first endpoint is within 1000 units of P
                %         %if pdist2(A,P)<1000
                %         % N is the vector pointing from A to the other end of the
                %         % link
                %         N = nodes(links(k,2)).coords - A;
                %         normN = norm(N); % the length of N
                %         % make the unit vector in the direction of N
                %         unitN = N/normN;
                %         % see equation at https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
                %         % The equation of a line can be given in vector form: X = A + t*unitN
                %         % in our case t is given by dot(P-A,unitN). But we don't want t to be
                %         % negative or to be bigger than normN or else it won't be within the line segment
                %         closestPointOnLine = A + min(normN,max(0,dot(P-A,unitN)))*unitN;
                %         distance = pdist2(closestPointOnLine, P);
                %         if distance < minDist
                %             closestPoint = closestPointOnLine;
                %             minDist = distance;
                %             closestLink = k;
                %         end
                %         %end
                %     end
            end
            % toc
        end
        
        function result = getLinkIndices(obj, nodeArray)
            % finds the indices for the links between nodes
            % nodeArray is an n x 2 array where each row is a pair of node
            % indices
            % returns an n x 1 array where the i'th element is the index of
            % the link between the nodes in the i'th row of nodeArray
            % the i'th element will be 0 if there is no such link (and will
            % print out a warning)
            % e.g.
            % nodeArray = [1 2; 3 4; 1 3];
            % getLinkIndices(nodeArray) returns the indices for the links
            % between nodes 1&2, 3&4, 1&3
            n = size(nodeArray,1);
            result = zeros(n,1);
            for i=1:n
                node1 = nodeArray(i,1);
                node2 = nodeArray(i,2);
                [~, linkIndex1] = ismember([node1 node2], obj.links(:,[1 2]), 'rows');
                [~, linkIndex2] = ismember([node2 node1], obj.links(:,[1 2]), 'rows');
                linkIndex = max(linkIndex1, linkIndex2);
                if ~linkIndex
                    disp(['getLinkIndices warning: returning 0 in element ' num2str(i)]);
                end
                result(i) = linkIndex;
            end
        end
        
        function obj = getSynSetGraph(obj, synSetIndex)
            % Creates the corresponding synGraph for the synSet defined by
            % synSetIndex (i.e. 1st, 2nd, 3rd...)
            %
            % make an undirected graph with all the synapses in synSet,
            % storing the distances between each connected pair of synapses
            % along the sekelton
            %
            % every synapse will have a list of the closest other synapses on either
            % side. If on a straight line alone, this will be 2 other synapse. But if
            % when you travel in one direction you hit a branch point before hitting
            % any synapses, you might have 2 or more other synapses in that direction.
            % 
            % the synGraph is a numSynapses x 1 cell array, in which the cell array
            % has an mx3 matrix (m is number of neighbouring synapses)
            % 1st column: index of the neighbouring synapse
            % 2nd column: distances to those neighbouring synapses
            % 3rd column: electrotonic distances to those neighbouring synapses
            % (2*distance*(sqrt(radius1)-sqrt(radius2))/(radius1-radius2) in each segment)
            % this is an "adjacency list" representation of an undirected graph
            % Then we can use Dijkstra's algorithm to calculate the distance from this
            % synapse to every other synapse
            
            disp('Finding distance between synapses');
            datestr(clock)
            
            s = obj.synSet(synSetIndex);
            
            numSyn = length(s.synLinks);
            
            % both real distance
            % and electrotonic distance
            synSetGraph = cell(numSyn,1);
            times = zeros(numSyn,1);
            WaitMessage = parfor_wait(numSyn, 'Waitbar', true, 'ReportInterval', 100);
            parfor i=1:numSyn
                tic
                synLoc = s.synLocs(i,:);
                % get which link the synapse is on
                synLink = s.synLinks(i);
                synSetGraph{i} = obj.synapseNeighbors(synSetIndex, synLoc, synLink, i);
                t = toc;
                fprintf("synapse %d took %f seconds\n", i, t);
                times(i) = t;
                WaitMessage.Send;
            end
            WaitMessage.Destroy;
            obj.synGraphs{synSetIndex} = synSetGraph;
            save(strcat(obj.synSetNames{synSetIndex},'-synDistGraph',datestr(now,'yyyymmddTHHMMSS'),'.mat'),'synSetGraph');
            save(strcat(obj.synSetNames{synSetIndex},'-synDistGraph',datestr(now,'yyyymmddTHHMMSS'),'-times.mat'),'times');
            
            datestr(clock)
        end
        
        function result = getNeighborsOnOtherGraph(obj, synSetIndices)
            % for each node on in graph synSetIndices(1), find the nearest
            % neighbors (and relevant distances) on graph synSetIndices(2)
            disp('Finding distance between synapses');
            datestr(clock)
            
            is = synSetIndices;
            % source = obj.synSet(is(1)); % commented out: no use sending
            % synSet and obj to the parfor workers
            
            numSource = length(obj.synSet(is(1)).synLinks);
            % both real distance
            % and electrotonic distance
            result = cell(numSource,1);
            times = zeros(numSource,1);
            
            WaitMessage = parfor_wait(numSource, 'Waitbar', true, 'ReportInterval', 100);
            parfor i=1:numSource
                tic
                synLoc = obj.synSet(is(1)).synLocs(i,:);
                % get which link the synapse is on
                synLink = obj.synSet(is(1)).synLinks(i);
                result{i} = obj.synapseNeighbors(is(2), synLoc, synLink, -1);
                t = toc;
                fprintf("synapse %d took %f seconds\n", i, t);
                times(i) = t;
                WaitMessage.Send;
                
            end
            save(strcat(obj.synSetNames{is(1)},'NeighborsOn',obj.synSetNames{is(2)},'-graph',datestr(now,'yyyymmddTHHMMSS'),'.mat'),'result');
            save(strcat(obj.synSetNames{is(1)},'NeighborsOn',obj.synSetNames{is(2)},'-graph',datestr(now,'yyyymmddTHHMMSS'),'-times.mat'),'times');
            
            datestr(clock)
            WaitMessage.Destroy
            
        end
        
        function result = synapseNeighbors(obj, synSetIndex, synLoc, synLink, synIndex)
            % return an n x 3 array
            % % 1st column: index of the neighbouring synapse
            % % 2nd column: distances to those neighbouring synapses
            % % 3rd column: electrotonic distances to those neighbouring synapses
            % % (2*distance*(sqrt(radius1)-sqrt(radius2))/(radius1-radius2) in each segment)
            
            % get the synapse's location on the skeleton
            s = obj.synSet(synSetIndex);
            
            result = [];

            % allow for rounding errors in comparing distances
            % Note: setting the distanceTolerance is important for avoiding
            % errors in finding synapse neighbors. 
            distanceTolerance = 0.001;
            
            % for both ends of the current link
            for j=1:2
                % first, ask are there any other synapses on the current link? if so,
                % then pick the nearest ones in each direction!
                % this needs to go inside the for loop because reused
                % the variable name synapseOnThisLink below - possible
                % mistake
                synapsesOnThisLink = find(s.synLinks==synLink);

                % have to add any synapses that might be on the end of the
                % link! eg on nextSkeletonNode
                % get the skeleton node in this direction
                nextSkeletonNode = obj.skel{obj.links(synLink,j),{'x','y','z'}};
                % search through synLocs for the 3d loc of nextSkeletonNode                
                [match,synapseOnNextNode] = ismember(nextSkeletonNode, s.synLocs, 'rows' );
                if match
                    synapsesOnThisLink = [synapsesOnThisLink; synapseOnNextNode];
                    synapsesOnThisLink = unique(synapsesOnThisLink);
                end

                
                % exclude the current synapse
                synapsesOnThisLink(synapsesOnThisLink==synIndex)=[];
                % exclude synapses that have the same location as synLoc
                % this is essential because this function stops once you've
                % found the closest neighbor. Thus if there are two
                % synapses in the same location, they will only see each
                % other as the closest neighbor. Thus they will be
                % disconnected from other nodes - a long chain of connected
                % nodes will stop with them.
                % But it creates another problem if you have two synapses
                % on top of each other: one of the pair will always get
                % ignored (no one will connect to it). So the true solution
                % not to have any duplicate synapse locations.
                % This is why it's crucial to combine all synapses that map
                % to the same location on the skeleton into a single
                % "synapse location"
                distancesOnThisLink = pdist2(synLoc,s.synLocs(synapsesOnThisLink,:));
                synapsesOnThisLink(distancesOnThisLink==0) = [];
                synLocToN = norm(nextSkeletonNode - synLoc);
                % if there were any synapses on this link, see if they lie on the line between
                % n and synLoc
                if ~isempty(synapsesOnThisLink)
                    distancesOnThisLink = pdist2(synLoc,s.synLocs(synapsesOnThisLink,:));
                    distancesToN = pdist2(nextSkeletonNode,s.synLocs(synapsesOnThisLink,:));
                    % use the triangle inequality
                    synapsesInThisDirection = find(distancesOnThisLink + distancesToN - synLocToN < distanceTolerance);
                    % note synapsesInThisDirection indexes within distancesOnThisLink
                    if ~isempty(synapsesInThisDirection)
                        % find the closest one
                        [minDist,index] = min(distancesOnThisLink(synapsesInThisDirection));
                        closestSynapseIndex = synapsesOnThisLink(synapsesInThisDirection(index));
                        % calculate the electrotonic distance to the
                        % closest synapse
                        % this calculation treats the segment as a
                        % truncated cone; see Methods of Amin et al 2020
                        % radius2 is the radius of the node in the
                        % 'forward' direction (i.e. the jth node of the
                        % link)
                        radius2 = obj.skel.radius(obj.links(synLink,j));
                        % radius1 is the radius in the 'backward' direction
                        % ie 1st node if j==2, 2nd node if j==1
                        radius1 = obj.skel.radius(obj.links(synLink,~(j-1)+1));
                        if (radius1==radius2)
                            minED = minDist / sqrt(radius1);
                        else
                            node1loc = obj.nodes(obj.links(synLink,~(j-1)+1)).coords;
                            linkDist = obj.linkDists(synLink,1);
                            x1 = pdist2(synLoc, node1loc);
                            x2 = pdist2(s.synLocs(closestSynapseIndex,:), node1loc);
                            radiusSlope = (radius2 - radius1)/linkDist;
                            minED = 2/radiusSlope*(sqrt(radiusSlope*x2 + radius1) - ...
                                sqrt(radiusSlope*x1 + radius1));
                        end
                        if minED<-distanceTolerance % allow for rounding errors
                            disp(strcat('node was ', num2str(obj.links(synLink,~(j-1)+1))));
                            disp(strcat('radius1 was ', num2str(radius1)));
                            disp(strcat('radius2 was ', num2str(radius2)));
                            disp(strcat('link was ', num2str(synLink)));
                            disp(strcat('synapse number was ', num2str(synIndex)));
                            error('electrotonic distance < 0');
                        end
                        if nnz(imag(minED))
                            disp(strcat('node was ', num2str(obj.links(synLink,~(j-1)+1))));
                            disp(strcat('radius1 was ', num2str(radius1)));
                            disp(strcat('radius2 was ', num2str(radius2)));
                            disp(strcat('link was ', num2str(synLink)));
                            disp(strcat('synapse number was ', num2str(synIndex)));
                            error('complex electrotonic distance');
                        end

                        % set up the next undirected node for the synapse graph for
                        % the index-th synapse
                        % that is, add the closest synapse to the list of this
                        % synapse's neighbours
                        result = [result; ...
                            [closestSynapseIndex, minDist, minED]];
                        continue
                    end
                end
                % if we have reached this point, we didn't find another synapse in
                % this direction, so we have to continue past the next skeleton
                % node
                
                % start a depth-first search for new synapses
                
                % push the index of the next node onto the stack
                nodeStack = obj.links(synLink,j);
                distanceStack = synLocToN;
                % calculate the electrotonic distance to the
                % next node, as above
                % radius2 is the radius of the node in the
                % 'forward' direction (i.e. the jth node of the
                % link)
                radius2 = obj.skel.radius(obj.links(synLink,j));
                % radius1 is the radius in the 'backward' direction
                % ie 1st node if j==2, 2nd node if j==1
                radius1 = obj.skel.radius(obj.links(synLink,~(j-1)+1));
                if (radius1==radius2)
                    EDtopush = synLocToN / sqrt(radius1);
                else
                    node1loc = obj.nodes(obj.links(synLink,~(j-1)+1)).coords;
                    linkDist = obj.linkDists(synLink,1);
                    x1 = pdist2(synLoc, node1loc);
                    % x2 = linkDist;
                    radiusSlope = (radius2 - radius1)/linkDist;
                    EDtopush = 2/radiusSlope*(sqrt(radius2) - ...
                        sqrt(radiusSlope*x1 + radius1));
                end
                if EDtopush<-distanceTolerance % allow for rounding errors
                    disp(strcat('node was ', num2str(obj.links(synLink,~(j-1)+1))));
                    disp(strcat('radius1 was ', num2str(radius1)));
                    disp(strcat('radius2 was ', num2str(radius2)));
                    disp(strcat('link was ', num2str(synLink)));
                    disp(strcat('synapse number was ', num2str(synIndex)));
                    error('electrotonic distance < 0');
                end
                if nnz(imag(EDtopush))
                    disp(strcat('node was ', num2str(obj.links(synLink,~(j-1)+1))));
                    disp(strcat('radius1 was ', num2str(radius1)));
                    disp(strcat('radius2 was ', num2str(radius2)));
                    disp(strcat('link was ', num2str(synLink)));
                    disp(strcat('synapse number was ', num2str(synIndex)));
                    error('complex electrotonic distance');                    
                end
                EDStack = EDtopush;  %electrotonic distance
                linksVisitedOnThisSearch = synLink;
                
                while ~isempty(nodeStack)
                    % nIndex <- nodeStack.pop()
                    nIndex = nodeStack(1);
                    nodeStack = nodeStack(2:end);
                    dist = distanceStack(1);
                    distanceStack = distanceStack(2:end);
                    ED = EDStack(1);
                    EDStack = EDStack(2:end);
                    
                    % get the links for this node
                    % in the first iteration of the while loop, the node will be links(synLink,j)
                    
                    neighbourNodesFromN = obj.nodes(nIndex).links;
                    linksFromN = obj.getLinkIndices([repmat(nIndex, length(neighbourNodesFromN), 1), neighbourNodesFromN(:)]);
                    % get rid of links that have already been visited
                    linksFromN = setdiff(linksFromN, linksVisitedOnThisSearch);
                    
                    % if this node has no unvisited links, then we are done with this node!
                    if isempty(linksFromN)
                        continue;
                    end
                    
                    % pick the first unvisited link
                    thisLink = linksFromN(1);
                    % mark it as visited
                    linksVisitedOnThisSearch = [linksVisitedOnThisSearch; thisLink];
                    
                    % if there is a synapse, we are done! Then we need to
                    % add up all the distances that we traversed to get
                    % here
                    synapsesOnThisLink = find(s.synLinks==thisLink);
                    % have to add any synapses that might be on the end of the
                    % link! eg on nextSkeletonNode
                    % get the skeleton node in this direction
                    nextSkeletonNode = obj.skel{obj.links(thisLink,j),{'x','y','z'}};
                    % search through synLocs for the 3d loc of nextSkeletonNode
                    [match,synapseOnNextNode] = ismember(nextSkeletonNode, s.synLocs, 'rows' );
                    if match
                        synapsesOnThisLink = [synapsesOnThisLink; synapseOnNextNode];
                        synapsesOnThisLink = unique(synapsesOnThisLink);
                    end
                    % get their distances from node nIndex
                    if ~isempty(synapsesOnThisLink)
                        distancesOnThisLink = pdist2(obj.nodes(nIndex).coords,s.synLocs(synapsesOnThisLink,:));
                        [minDist,index] = min(distancesOnThisLink);
                        closestSynapseIndex = synapsesOnThisLink(index);
                        % calculate the electrotonic distance from node nIndex to the
                        % closest synapse, as above
                        thisLinksNodes = obj.links(thisLink,:);
                        % radius2 is the radius of the node in the
                        % 'forward' direction (i.e. not the node nIndex)
                        radius2 = obj.skel.radius(thisLinksNodes(thisLinksNodes~=nIndex));
                        % radius1 is the radius of the node where you are
                        % right now
                        radius1 = obj.skel.radius(nIndex);
                        if (radius1==radius2)
                            minED = minDist / sqrt(radius1);
                        else
                            node1loc = obj.nodes(nIndex).coords;
                            linkDist = obj.linkDists(thisLink,1);
                            % x1 = 0
                            x2 = pdist2(s.synLocs(closestSynapseIndex,:), node1loc);
                            radiusSlope = (radius2 - radius1)/linkDist;
                            minED = 2/radiusSlope*(sqrt(radiusSlope*x2 + radius1) - sqrt(radius1));
                        end
                        if minED<-distanceTolerance % allow for rounding errors
                            disp(strcat('node was ', num2str(nIndex)));
                            disp(strcat('radius1 was ', num2str(radius1)));
                            disp(strcat('radius2 was ', num2str(radius2)));
                            disp(strcat('link was ', num2str(thisLink)));
                            disp(strcat('synapse number was ', num2str(synIndex)));
                            error('electrotonic distance < 0');
                        end
                        if nnz(imag(minED))
                            disp(strcat('node was ', num2str(nIndex)));
                            disp(strcat('radius1 was ', num2str(radius1)));
                            disp(strcat('radius2 was ', num2str(radius2)));
                            disp(strcat('link was ', num2str(thisLink)));
                            disp(strcat('synapse number was ', num2str(synIndex)));
                            error('complex electrotonic distance');
                        end
                        result = [result; ...
                            [closestSynapseIndex, dist + minDist, ED + minED]];
                        
                    end
                    
                    if isempty(synapsesOnThisLink)
                        % need to push the next node
                        % but also push the current one in case it
                        % still has other links to explore
                        theseNodes = obj.links(thisLink,1:2);
                        nextNode = setdiff(theseNodes, nIndex);
                        
                        % Need to push different distances on the stack
                        % for each node: for the next node, add on the
                        % current link. But for the current link, do
                        % not add on the current link.
                        if length(linksFromN)>1
                            % push the current node on the bottom of the stack so we continue a depth first search
                            nodesToPush = [nextNode; nIndex];
                            distToPush = [dist + obj.linkDists(thisLink,1); dist];
                            EDToPush = [ED + obj.linkDists(thisLink,2); ED];
                        else
                            nodesToPush = nextNode;
                            distToPush = dist + obj.linkDists(thisLink,1);
                            EDToPush = ED + obj.linkDists(thisLink,2);
                        end
                        % These next options are if you did find a synapse
                        % on the link
                    elseif length(linksFromN)>1
                        % If there are more links from this node, you
                        % are not done with this node yet! you
                        % have to push the node back onto the stack
                        nodesToPush = nIndex;
                        distToPush = dist;
                        EDToPush = ED;
                        % If there are no other links, don't push anything
                        % on the stack - the next runs on the while loop
                        % will deal with the remaining nodes
                    else
                        nodesToPush = [];
                        distToPush = [];
                        EDToPush = [];
                    end
                    
                    distanceStack = [distToPush; distanceStack];
                    EDStack = [EDToPush; EDStack];
                    nodeStack = [nodesToPush; nodeStack];
                end
            end % for j=1:2
        end

        function [distances] = getAllSynDistWithinSynSet(obj, synSetIndex, maxDist, whichDist)
            % Use Dijkstra's algorithm to find the distance from every
            % synapse to every other synapse in a synSet (given by
            % synSetIndex)
            % maxDist is to set an upper limit on distances - if you know
            % the skeleton only has a total length of X, it's impossible
            % for the distances between synapse to be greater than X
            % whichDist should be 2 for regular distance, 3 for electrotonic distance 
            
            graph = obj.synGraphs{synSetIndex};
            numNodes = length(graph);
            
            % make square empty arrays
            distances = Inf(numNodes, numNodes, 'single');
            %prevNodes = single(numNodes, numNodes, 'single');
            
            WaitMessage = parfor_wait(numNodes, 'Waitbar', true, 'ReportInterval', 10);
            
            parfor j=1:numNodes
                % we are going to loop over ROWS
                theseDistances = Inf(numNodes, 1, 'single');
                % Implements the pseudocode from the Wikipedia page about
                % Dijkstra's algorithm
                theseDistances(j) = 0; % distance of j to itself is 0
                Q = 1:numNodes; % you don't have to measure distances from 1-j because
                % they've already been measured on iterations before j!
                while ~isempty(Q)
                    % u is the index of the node with the minimum distance
                    % but only search for the distances of nodes that are still
                    % in Q
                    [minDist, minDistIndex] = min(theseDistances(Q));
                    if (minDist > maxDist)
                        % if even the closest distance is already over the
                        % limit, there is no way we will find any other nodes
                        % with a distance smaller than maxDist. So give up!
                        break;
                    end
                    u = Q(minDistIndex);
                    % remove u from Q
                    Q = Q(Q~=u); % ie keep only values in Q that are not u
                    
                    % for each neighbor of u
                    for i=1:size(graph{u},1)
                        % v is the index of the neighbor
                        v = graph{u}(i,1);
                        % alt is the distance of node u, plus the distance
                        % between u and the neighbor
                        
                        alt = theseDistances(u) + graph{u}(i,whichDist);
                        if (alt < theseDistances(v)) % if we've found a shorter distance for v
                            theseDistances(v) = alt;
                            %prevNodes(v) = u;
                        end
                    end
                end
                distances(j,:) = theseDistances;
                WaitMessage.Send;

            end
            
            WaitMessage.Destroy
        end

        function [] = drawSkeleton(obj)
            % This is for real neurons but beware of very large skeletons
            
            % ignore links that have 0 values
            linkIndices = 1:size(obj.links,1);
            nonzeroLinkIndices = linkIndices(all(obj.links,2));
            if length(nonzeroLinkIndices)>10000
                resp = input('This is a huge skeleton (>10000 nodes). Are you sure you want to continue? press n to exit');
                if (resp=='n')
                    return;
                end
            end
            
            figure
            hold on;
            for j=nonzeroLinkIndices
                if mod(j,1000)==0
                    j
                end
                coords(1,:) = obj.nodes(obj.links(j,1)).coords; % 1x3
                coords(2,:) = obj.nodes(obj.links(j,2)).coords; % 1x3
                % coords is now a 2x3 array
                if ~isempty(obj.linksInROI) && obj.linksInROI(j)
                    line(coords(:,2)',coords(:,1)',-coords(:,3)','Color','r','LineWidth',1);
                else
                    line(coords(:,2)',coords(:,1)',-coords(:,3)','Color','k','LineWidth',1);
                end
            end
%             nodeIndices = 1:length(obj.nodes);
%             for i=nodeIndices
%                 x = obj.nodes(i).coords(2);
%                 y = obj.nodes(i).coords(1);
%                 z = -obj.nodes(i).coords(3);
%                 plot3(x,y,z,'o','Color','k'); % need to change this when having multiple synSets!
%             end
%             cmap = colormap;
            colorList = ['m','c','g','r','b','k','y'];
            for i=1:length(obj.synSet)
                color = colorList(i); %cmap(i*32,:);
                for j=1:length(obj.synSet(i).synLinks)
                    x = obj.synSet(i).synLocs(j,2);
                    y = obj.synSet(i).synLocs(j,1);
                    z = -obj.synSet(i).synLocs(j,3);
                    plot3(x,y,z,'o','MarkerSize',3,'MarkerEdgeColor',color,'MarkerFaceColor',color); % need to change this when having multiple synSets!
                    %text(x,y,z,strcat('...',num2str(j)));
                end
                % number of synapses
%                 for j=1:length(obj.synGraphs{i})
%                     % number of neighbors
%                     for k=1:size(obj.synGraphs{i}{j},1)
%                         coords(1,:) = obj.synSet(i).synLocs(j,:);
%                         coords(2,:) = obj.synSet(i).synLocs(obj.synGraphs{i}{j}(k,1),:);
%                         line(coords(:,1)',coords(:,2)',coords(:,3)','Color',color,'LineWidth',1);
%                     end
%                 end
            end
            
%             linksPerNode = obj.linksPerNode();
%             nodesWithBranches = find(linksPerNode>2);
%             
%             for i=1:length(nodesWithBranches)
%                 loc = obj.nodes(nodesWithBranches(i)).coords;
%                 plot3(loc(2),loc(1),-loc(3),'o','Color','k');
%             end
            axis equal, axis tight
        end
        
        function obj = labelInROI(obj, ROIvertices, ROIfaces)
            % ROI is a mesh
            % this can't be the final code - need to make it general to
            % allow multiple ROIs to be labeled?
            nodeLocs = table2array(obj.skel(:,{'x','y','z'}));
            obj.nodesInROI = intriangulation(ROIvertices, ROIfaces, nodeLocs);
            
            obj.linksInROI = false(size(obj.links,1),1);
            % ignore links that have 0 values
            linkIndices = 1:size(obj.links,1);
            nonzeroLinkIndices = linkIndices(all(obj.links,2));
            nodes1 = obj.links(nonzeroLinkIndices,1);
            nodes2 = obj.links(nonzeroLinkIndices,2);
            obj.linksInROI(nonzeroLinkIndices) = obj.nodesInROI(nodes1) & obj.nodesInROI(nodes2);
        end

        
        function [] = drawTestSkeleton(obj)
            % This is not for real neurons, but for debugging - use it on
            % the toy skeletons to make sure that synapses are finding the
            % right neighbors, being mapped correctly onto the skeleton,
            % etcv
            
            % ignore links that have 0 values
            linkIndices = 1:size(obj.links,1);
            nonzeroLinkIndices = linkIndices(all(obj.links,2));
            if length(nonzeroLinkIndices)>10000
                resp = input('This is a huge skeleton (>10000 nodes). Are you sure you want to continue? press n to exit');
                if (resp=='n')
                    return;
                end
            end
            
            figure
            hold on;
            for j=nonzeroLinkIndices
                if mod(j,1000)==0
                    j
                end
                coords(1,:) = obj.nodes(obj.links(j,1)).coords; % 1x3
                coords(2,:) = obj.nodes(obj.links(j,2)).coords; % 1x3
                % coords is now a 2x3 array
                line(coords(:,1)',coords(:,2)',coords(:,3)','Color','k','LineWidth',1);
            end
            nodeIndices = 1:length(obj.nodes);
            for i=nodeIndices
                x = obj.nodes(i).coords(1);
                y = obj.nodes(i).coords(2);
                z = obj.nodes(i).coords(3);
                plot3(x,y,z,'o','Color','k'); % need to change this when having multiple synSets!
                text(x,y,z,strcat('...',num2str(i)));
            end
            cmap = colormap;
            for i=1:length(obj.synSet)
                color = cmap(i*32,:);
                for j=1:length(obj.synSet(i).synLinks)
                        x = obj.synSet(i).synLocs(j,1);
                    y = obj.synSet(i).synLocs(j,2);
                    z = obj.synSet(i).synLocs(j,3);
                    plot3(x,y,z,'o','MarkerSize',5,'MarkerEdgeColor',color,'MarkerFaceColor',color); % need to change this when having multiple synSets!
                    text(x,y,z,strcat('...',num2str(j)));
                end
                % number of synapses
                for j=1:length(obj.synGraphs{i})
                    % number of neighbors
                    for k=1:size(obj.synGraphs{i}{j},1)
                        coords(1,:) = obj.synSet(i).synLocs(j,:);
                        coords(2,:) = obj.synSet(i).synLocs(obj.synGraphs{i}{j}(k,1),:);
                        line(coords(:,1)',coords(:,2)',coords(:,3)','Color',color,'LineWidth',1);
                    end
                end
            end
            axis equal, axis tight
        end
        
        
        function [] = drawSparseSynapses(obj, numToDraw, j, varargin)
            % use this to only draw a subsample all the synapses on a
            % neurons to avoid overwhelming Matlab
            figure; hold on;
            % use the commented-out code to color each synapse according to
            % which ROI it's tagged with
            if nargin>3
                inputColor=varargin{1};
            else
                inputColor='k';
            end
            cmap = colormap;
            synsToUse = randsample(length(obj.synSet(j).synLinks), numToDraw);
            for i=1:length(synsToUse)
                x = .008*obj.synSet(j).synLocs(synsToUse(i),2);
                y = .008*obj.synSet(j).synLocs(synsToUse(i),1);
                z = -.008*obj.synSet(j).synLocs(synsToUse(i),3);
                scatter3(x,y,z,'o','MarkerEdgeColor',inputColor,'MarkerEdgeAlpha',0.3, 'LineWidth',0.5); % 'Color',color); %'MarkerSize',5, 
            end
            ROIstoHighlight = {'','SLP(R)','SCL(R)','ICL(R)','PLP(R)'}; %,'CRE(R)','SIP(R)'
            for i=1:length(ROIstoHighlight)
                ROIindex = find(strcmp(ROIstoHighlight(i), obj.synSet(j).synROIskey));
                % Test the calyxHACK
                synsInThisROI = (obj.synSet(j).synROIs==ROIindex)&(obj.synSet(j).synLocs(:,2)<(160/.008));
                x = .008*obj.synSet(j).synLocs(synsInThisROI,2);
                y = .008*obj.synSet(j).synLocs(synsInThisROI,1);
                z = -.008*obj.synSet(j).synLocs(synsInThisROI,3);
                color = cmap(floor(i*63/length(ROIstoHighlight))+1,:);
                scatter3(x,y,z,'o','MarkerEdgeColor',color);
                scatter3(180,180,-40-i*10,'o','MarkerEdgeColor',color);
                text(180,180,-40-i*10,strcat('...',ROIstoHighlight(i)));
            end
            axis equal, axis tight;
            view(66,15)
        end
        
        function [] = drawSparseSynapsesAPL(obj, numToDraw, j)
            % This is for Fig 8C - re-indexing the Voronoi segments of the
            % hemibrain APL to start from the calyx and go toward the lobes
            figure; hold on;
            cmap = colormap;
            maxSegment = max(obj.synSet(j).voronoiSegments);
            distFromCalyx = [19	20	18	20	21	17	21	22	16	22	23	15	23	24	14	24	25	13	25	26	12	26	27	11	28	10	9	8	7	6	5	4	3	2	1];
            synsToUse = find(ismember(obj.synSet(j).synROIs, [3 7 13:20]) & (obj.synSet(j).synLocs(:,1) > 10000) & (obj.synSet(j).synLocs(:,1) < 28000) & obj.synSet(j).voronoiSegments);
%             for i=randsample(length(obj.synSet(j).synLinks), numToDraw)
            synsToUse = randsample(synsToUse, numToDraw);
            for i=1:length(synsToUse)
                x = .008*obj.synSet(j).synLocs(synsToUse(i),2);
                y = .008*obj.synSet(j).synLocs(synsToUse(i),1);
                z = -.008*obj.synSet(j).synLocs(synsToUse(i),3);
                try
                    color = cmap(floor(distFromCalyx(obj.synSet(j).voronoiSegments(synsToUse(i)))*63/max(distFromCalyx))+1,:);
                catch
%                     distFromCalyx(obj.synSet(j).voronoiSegments(i))
                    color = [0 0 0];
                end
                plot3(x,y,z,'o','MarkerFaceColor',color,'MarkerEdgeColor','k','MarkerSize',5);
            end
            axis equal, axis tight;
            view(66,15);
        end
        
        function [] = drawKCSkeleton(obj)
            
            % ignore links that have 0 values
            linkIndices = 1:size(obj.links,1);
            nonzeroLinkIndices = linkIndices(all(obj.links,2));
            if length(nonzeroLinkIndices)>10000
                resp = input('This is a huge skeleton (>10000 nodes). Are you sure you want to continue? press n to exit');
                if (resp=='n')
                    return;
                end
            end
            
            figure
            hold on;
            % draw skeleton links
            for j=nonzeroLinkIndices
                if mod(j,1000)==0
                    j
                end
                coords(1,:) = obj.nodes(obj.links(j,1)).coords*0.008; % 1x3
                coords(2,:) = obj.nodes(obj.links(j,2)).coords*0.008; % 1x3
                line(coords(:,2)',coords(:,1)',-coords(:,3)','Color','k','LineWidth',2);
            end
            [maxRadius, biggestNodeIndex] = max(obj.skel.radius);
            if maxRadius>100
                somaCoords = 0.008*[1 1 -1].*obj.nodes(biggestNodeIndex).coords;
                x = somaCoords(2);
                y = somaCoords(1);
                z = somaCoords(3);
                plot3(x,y,z,'o','MarkerSize',15,'MarkerEdgeColor','k','MarkerFaceColor','k');
            end
            % draw APL synapses
            APLtoKC = find(obj.synSet(2).uniquePartners == 425790257);
            APLtoKCsyns = find(obj.synSet(2).uniquePartnerIndices == APLtoKC);
            % The commented out code allows you to color-code synapses by
            % ROI
            % This reveals that many calyx synapses are annotated as
            % belonging to other ROIs besides CA(R)
%             calyxIndex = find(strcmp('CA(R)', obj.synSet(2).synROIskey));
%             cmap = colormap;
%             ROIsdisplayed=[];
            for j=1:length(APLtoKCsyns)
                thisSyn = obj.synSet(2).synLocsOfOrigSyn(APLtoKCsyns(j));
%                 thisROI = obj.synSet(2).synROIs(thisSyn);
%                 if ~ismember(thisROI, ROIsdisplayed)
%                     ROIsdisplayed = [ROIsdisplayed, thisROI];
%                 end    
%                 color = cmap(floor(obj.synSet(2).synROIs(thisSyn)*63/length(obj.synSet(2).synROIskey))+1,:);
                x = 0.008*obj.synSet(2).synLocs(thisSyn,2);
                y = 0.008*obj.synSet(2).synLocs(thisSyn,1);
                z = -0.008*obj.synSet(2).synLocs(thisSyn,3);
%                 if ~sum(strcmp(obj.synSet(2).synROIskey(obj.synSet(2).synROIs(thisSyn)), {'CA(R)','PED(R)','a''L(R)','aL(R)','b''L(R)','bL(R)','gL(R)'}))
%                     obj.synSet(2).synROIskey(obj.synSet(2).synROIs(thisSyn))
%                     text(x,y,z,strcat('...', obj.synSet(2).synROIskey(obj.synSet(2).synROIs(thisSyn))));
%                 end
%                 if obj.synSet(2).synROIs(obj.synSet(2).synLocsOfOrigSyn(APLtoKCsyns(j)))==calyxIndex
%                     plot3(x,y,z,'o','MarkerSize',6,'MarkerEdgeColor','g','MarkerFaceColor','g');
%                 else
                plot3(x,y,z,'o','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','r'); %replace 'r' with color to color-code synapses
%                 end
            end
            KCtoAPL = find(obj.synSet(1).uniquePartners == 425790257);
            KCtoAPLsyns = find(obj.synSet(1).uniquePartnerIndices == KCtoAPL);
            for j=1:length(KCtoAPLsyns)
                x = 0.008*obj.synSet(1).synLocs(obj.synSet(1).synLocsOfOrigSyn(KCtoAPLsyns(j)),2);
                y = 0.008*obj.synSet(1).synLocs(obj.synSet(1).synLocsOfOrigSyn(KCtoAPLsyns(j)),1);
                z = -0.008*obj.synSet(1).synLocs(obj.synSet(1).synLocsOfOrigSyn(KCtoAPLsyns(j)),3);
                plot3(x,y,z,'o','MarkerSize',6,'MarkerEdgeColor','b','MarkerFaceColor','b');
            end
            axis equal, axis tight
            view(66,15);
            xlim([70  290])
            ylim([57  220])
            zlim([-182.0300  -34.5])
%             for i=1:length(ROIsdisplayed)
%                 color = cmap(floor(ROIsdisplayed(i)*63/length(obj.synSet(2).synROIskey))+1,:);
%                 plot3(100,100,-182+i*10,'o','MarkerSize',6,'MarkerEdgeColor',color,'MarkerFaceColor',color);
%                 text(100,100,-182+i*10,strcat('...', obj.synSet(2).synROIskey(ROIsdisplayed(i))));
%             end
            saveas(gcf,strcat(num2str(obj.synSet(1).rawData.upstream_bodyId(1)),'.png'),'png');
        end
        
        function numPartners = neighborsPerSynLoc(obj,j)
            % in synSet j, how many neighboring synapses does each synapse
            % have?
            numPartners = zeros(length(obj.synSet(j).synLinks),1);
            for i=1:size(numPartners,1)
                numPartners(i) = size(obj.synGraphs{j}{i},1);
            end
        end
        function result = synapsesPerSynLoc(obj,j)
            % in synSet j, how many neighboring synapses does each synapse
            % have?
            result = zeros(length(obj.synSet(j).synLinks),1);
            for i=1:size(result,1)
                result(i) = length(obj.synSet(j).synGroupings{i});
            end
        end
        function result = linksPerNode(obj)
            % for each node, how many links are there?
            result = zeros(length(obj.nodes),1);
            for i=1:length(obj.nodes)
                result(i) = length(obj.nodes(i).links);
            end
        end
        
        % next 4 functions assme there is already a figure open
        function [] = drawLink(obj,j)
            % assuming there is already a figure open
            coords(1,:) = obj.nodes(obj.links(j,1)).coords; % 1x3
            coords(2,:) = obj.nodes(obj.links(j,2)).coords; % 1x3
            % coords is now a 2x3 array
            line(coords(:,1)',coords(:,2)',coords(:,3)','Color','k','LineWidth',1);
            
        end
        
        function [] = drawSynLoc(obj,i,j,color)
            loc = obj.synSet(i).synLocs(j,:);
            plot3(loc(1),loc(2),loc(3),'o','Color',color);
            text(loc(1),loc(2),loc(3),strcat('...',num2str(j)));
            
        end
        
        function [] = drawNode(obj,j)
            loc = obj.nodes(j).coords;
            plot3(loc(1),loc(2),loc(3),'o','Color','k');
            text(loc(1),loc(2),loc(3),strcat('...',num2str(j)));
        end
        
        function [] = drawSynLink(obj,i,j)
            obj.drawLink(obj.synSet(i).synLinks(j));
            
        end
        
        function [] = drawSparseBranchingNodes(obj, numToDraw, branches)
            linksPerNode = obj.linksPerNode();
            nodesWithBranches = find(linksPerNode>branches);
            
            toDraw = randsample(length(nodesWithBranches),numToDraw);
            figure
            hold on
            for i=1:length(toDraw)
                loc = obj.nodes(nodesWithBranches(toDraw(i))).coords;
                plot3(loc(2),loc(1),-loc(3),'o','Color','k');
            end
        end
        
        function obj = mapOrigSynsToSynLocs(obj)
            % This function should also be unnecessary if we keep "ia" from
            % getSynSets            
            for i=1:length(obj.synSet)
                obj.synSet(i).synLocsOfOrigSyn = zeros(length(obj.synSet(i).uniquePartnerIndices),1,'single');
                for j=1:length(obj.synSet(i).synGroupings)
                    obj.synSet(i).synLocsOfOrigSyn(obj.synSet(i).synGroupings{j}) = j;
                end
            end
        end
        
        function obj = scrambleUniquePartners(obj, synSetIndices)
            % WARNING: this will IRREVERSIBLY scramble the unique partner
            % IDs of all synapses. save a copy of your object before you
            % run this, or put the result into a DIFFERENT variable!
            is = synSetIndices;
            for i=1:length(is)
                numPartners = length(obj.synSet(is(i)).uniquePartnerIndices);
                obj.synSet(is(i)).uniquePartnerIndices = ...
                    obj.synSet(is(i)).uniquePartnerIndices(randperm(numPartners));
            end
            
            
        end
        
        function result = weightedSumBetweenSynSets(obj, decayFactor, distMatrix, synSetIndices, ROIs, varargin)
            % distMatrix is KCtoAPLdisttoAPLtoKC
            % synSetIndices: dim 1 is the synSet for dim 1 of distMatrix,
            % and likewise for dim 2
            % we want to apply ROIs only to the output synapses of APL
            % (APLtoKC synapses) this should be dim 2 of distMatrix (synSet
            % 1)
            % columns of output correspond to dim 2 of distMatrix
            
            if nargin>5
                calyxHACK=varargin{1};
            else
                calyxHACK=false;
            end

            
            is = synSetIndices;
            
            ROItags = cell(length(ROIs),1);
            for i=1:length(ROItags)
                [~,ROItags{i},~] = intersect(obj.synSet(is(2)).synROIskey,ROIs{i});
%                 ROItags(i) = find(strcmp(obj.synSet(is(2)).synROIskey,ROIs{i}));
            end
            ROItags
            ROItags{1}
            tic
            if decayFactor==Inf
                weighting = true(size(distMatrix));
            else
                weighting = exp(-distMatrix/decayFactor);
            end
            toc
            numUniquePartners = [length(obj.synSet(is(1)).uniquePartners), ...
                length(obj.synSet(is(2)).uniquePartners)];
            result = zeros(numUniquePartners(1), numUniquePartners(2), 'single'); %length(ROIs), 'single');
            WaitMessage = parfor_wait(numUniquePartners(1), 'Waitbar', true, 'ReportInterval', 1);
            
           % this will end up in the row dimension (dim 1) of result
            for i=1:numUniquePartners(2) % loop over synapses in dim 2 of distMatrix. that's APLtoKC synapses
                % all the synapses for the ith unique partner
                %                 synapses1 = (obj.synSet(is(1)).uniquePartnerIndices==i);
                % get the synLocs of these synapses
                tic
                for r=1:length(ROItags)
                    % only apply ROI filter to dim 2 of distMatrix
                    origSynLocs = (obj.synSet(is(2)).uniquePartnerIndices==i); % find original synapses with partner i
                    synLocs2 = obj.synSet(is(2)).synLocsOfOrigSyn(origSynLocs); % find the uniqueSynapses of these original synapses
                    if calyxHACK
                        % This is a terrible hack for isolating the calyx!
                        % necessary because not all calyx synapses are
                        % tagged CA(R)
                        % Many non-CA(R) synapses are actually in
                        % the calyx, but some are not (some are in the peduncle, 
                        % some are in the lobes. So to include the
                        % non-CA(R) synapses and exclue the non-calyx
                        % synapses, throw away all synapses whose x value
                        % is > 160 um
                        synLocs2 = synLocs2((obj.synSet(is(2)).synLocs(synLocs2,2)<(160/.008)));
                    end
                    % reason for indexing here and not in the nested for loop
                    % is that Matlab is slow at extracting rows in a very large
                    % matrix. So we extract a few columns here to send to the
                    % for loop.
                    lia = ismember(obj.synSet(is(2)).synROIs(synLocs2), ROItags{r}); % find only those uniqueSynapses that match the ROItags
                    fprintf('ROI %d makes up fraction %f of synapses\n', r, nnz(lia)/length(synLocs2));
                    rowWeighting = weighting(:, synLocs2(lia));
                    % this will end up in the column dimension (dim 2) of result
                    for j=1:numUniquePartners(1) % loop over synapses in dim 1 of distMatrix
                        % all the synapses for the jth unique partner
                        %                     synapses2 = (obj.synSet(is(2)).uniquePartnerIndices==j);
                        % get the synLocs of these synapses
                        synLocs1 = obj.synSet(is(1)).synLocsOfOrigSyn((obj.synSet(is(1)).uniquePartnerIndices==j));
                        allWeights = rowWeighting(synLocs1, :);
                        result(i,j,r) = sum(allWeights(:));
                    end
                end
                t=toc;
                fprintf('row %d took %f s for %d synapses\n',i,t,length(synLocs2));
                WaitMessage.Send;
            end
            WaitMessage.Destroy;
        end
        
        function obj = addEvenRandomChans(obj, n, rIndex)
            % add n synapses randomly and evenly distributed throughout the
            % skeleton's links
            % only on the main skeleton
            % create this as synSet # rIndex (r = random)
            if (rIndex==1)||(rIndex==2)
                error('these synSets already exist!');
            end
            
            linksInMainSkel = find(obj.labeledLinks==mode(obj.labeledLinks));
            totalLength = 0;
            for i=linksInMainSkel'
                totalLength = totalLength+pdist2(obj.nodes(obj.links(i,1)).coords, obj.nodes(obj.links(i,2)).coords);
                if mod(i,10000)==0
                    fprintf('measuring length %d\n', i);
                end
            end
            
            obj.synSet(rIndex).synLinks = [];
            obj.synSet(rIndex).synLocs = [];
            
            % do not make random synapses outside the main skeleton
            for i=linksInMainSkel'
                if mod(i,10000)==0
                    fprintf('placing synapses %d\n', i);
                end
                if sum(obj.links(i,:))
                    A = obj.nodes(obj.links(i,1)).coords;
                    B = obj.nodes(obj.links(i,2)).coords;
                    linkLength = pdist2(A, B);
                    P = n/totalLength*linkLength;
                    count = floor(P);
                    
                    if rand(1)<rem(P,1)
                        count = count+1;
                    end
                    if count
                        for j=1:count
                            loc = A + (B-A)*rand(1);
                            obj.synSet(rIndex).synLocs = [obj.synSet(rIndex).synLocs; loc];
                            obj.synSet(rIndex).synLinks = [obj.synSet(rIndex).synLinks; i];
                        end
                    end
                end
            end
            
            obj.synSetNames{rIndex} = 'random';
        end
        
        function total = getTotalLengthMainBranch(obj, ED)
            % Get the total length of the branch of the skeleton that has
            % the most links
            % ED = 1 to get electrotonic distance, 0 to get actual distance
            if isempty(obj.labeledLinks)
                obj = obj.labelNodes();
            end
            
            linksInMainSkel = (obj.labeledLinks==mode(obj.labeledLinks));
            
            total = sum(obj.linkDists(linksInMainSkel,ED+1));
            
            % below is obsolete code from treating links as cylinders
            % instead of truncated cones
%             total = 0;
%             if ED
%                 for i=linksInMainSkel'
%                     total = total+pdist2(obj.nodes(obj.links(i,1)).coords, obj.nodes(obj.links(i,2)).coords)/sqrt(obj.links(i,3));
%                     if mod(i,1000)==0
%                         disp(i);
%                     end
%                 end
%             else
%                 for i=linksInMainSkel'
%                     total = total+pdist2(obj.nodes(obj.links(i,1)).coords, obj.nodes(obj.links(i,2)).coords);
%                     if mod(i,1000)==0
%                         disp(i);
%                     end
%                 end
%             end
        end
        
        function result = getNumSynapsesPerPartner(obj, i)
            numPartners = length(obj.synSet(i).uniquePartners);
            result = zeros(numPartners,1);
            for j=1:numPartners
                result(j) = nnz(obj.synSet(i).uniquePartnerIndices==j);
            end
        end

        function result = diffRadiiPerNode(obj)
            result = nan(size(obj.skel,1),1);
            for i=1:length(result)
                if obj.skel.link(i)>0
                    thisRadius = obj.skel.radius(i);
                    linkedRadius = obj.skel.radius(obj.skel.link(i));
                    result(i) = (thisRadius - linkedRadius)/(thisRadius + linkedRadius);
                end
            end
        end
        
        function meanRadius = meanRadiusOverall(obj)
            nodes1 = obj.links(2:end,1);
            nodes2 = obj.links(2:end,2);
            distances = obj.linkDists(2:end,1);
            radii = mean([obj.skel.radius(nodes1), obj.skel.radius(nodes2)], 2);
            meanRadius = sum(radii.*distances)/sum(distances); %weighted mean
        end
        
        
    end %methods
    
    methods (Static)
        function distances = getDistBetweenSynSets(distMatrix, neighbors, ED)
            % distMatrix: distances between every synapse in one synSet
            % (say synSet 1)
            % (square matrix)
            % neighbors: a cell array where each cell holds the neighbors
            % in synSet 1 of each synapse in synSet 2, and the distances
            % (the results of getNeighborsOnOtherGraph)
            % ED: 0 for real distance, 1 for electrotonic distance
            
            numSource = size(neighbors,1);
            numTarget = size(distMatrix,1);
            
            distances = zeros(numSource,numTarget,'single');
            datetime(clock)
            WaitMessage = parfor_wait(numSource, 'Waitbar', true, 'ReportInterval', 100);
            
            for i=1:numSource
%                 tic
                
                if ~isempty(neighbors{i})
                    [minDists, minDistIndices] = min(distMatrix(neighbors{i}(:,1),:),[],1);
                    distances(i,:) = minDists + neighbors{i}(minDistIndices,2+ED)';
                end
                WaitMessage.Send;
%                 t=toc;
%                 fprintf("Row %d took %f s\n", i,t);
                
            end
            datetime(clock)
            WaitMessage.Destroy
        end
    end
    
    
    
end