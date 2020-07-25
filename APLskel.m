classdef APLskel < neurSkel
    properties
        linkVoronoiSegments
        nodeVoronoiSegments
    end
    
    methods
        
        function obj = APLskel()
            obj = obj@neurSkel('APLskelHealed.csv');
        end
 
        
        function obj = getKCSynapses(obj)
            obj = obj.getSynSets({'APLtoKC.csv', 'KCtoAPL.csv'}, [1 2]);
        end
        
        function obj = assignSynapsesToVoronoiSegments(obj)
            % assign synapses to voronoi segments according to which link
            % they are on
            for i=1:length(obj.synSet)
                obj.synSet(i).voronoiSegments = obj.linkVoronoiSegments(obj.synSet(i).synLinks);
            end
        end
        
        function [result, activityPerSyn] = stimVoronoi(obj, decayFactor, distMatrix, i, stim)
            % distMatrix is the distance between random synapses
            % stim is a numVoronoiSegments x 1 list of red dye values at
            % each voronoi segment
            % i is which synSet (for APL, this is 3)
            tic; weights = exp(-distMatrix/decayFactor); toc
            if size(distMatrix,1)~=length(obj.synSet(i).synLinks)
                error('dimensions don''t match');
            end
            % this is to use only the random points within the calyx,
            % peduncle or lobes, and to ignore the points that are in
            % processes going outside the right mushroom body
            synsToUse = find(ismember(obj.synSet(i).synROIs, [3 7 13:18]) & ...
                (obj.synSet(i).synLocs(:,1) > 10000) & (obj.synSet(i).synLocs(:,1) < 28000) ...
                & obj.synSet(i).voronoiSegments);
            
            numSyn = length(synsToUse);
            stimPerSyn = zeros(numSyn,1);
            tic
            for j=1:numSyn
                % each synLoc is stimulated according to the red dye at
                % that segment
                vseg = obj.synSet(i).voronoiSegments(synsToUse(j));
                if vseg
                    stimPerSyn(j) = stim(vseg);
                end
            end
            toc
            activityPerSyn = zeros(numSyn,1);
            tic
            for j=1:numSyn
                % each synLoc is activated by all other synLocs, weighted
                % by exponential decay with distance
                activityPerSyn(j) = sum(stimPerSyn.*weights(synsToUse,synsToUse(j)));
            end
            toc
            numVoronoiSegments = max(obj.synSet(i).voronoiSegments);
            result = zeros(numVoronoiSegments,1);
            tic
            
            for j=1:numVoronoiSegments
                result(j) = mean(activityPerSyn(obj.synSet(i).voronoiSegments(synsToUse)==j));
            end
            toc
            
        end
        
        function obj = labelRandomSynsByROI(obj)
            % the real synapses need to be labeled already
            % label the random synapses according to the nearest synapse in
            % synSet 2 (the KC-APL synapses)
            labels = zeros(length(obj.synSet(3).synLinks), 1);
            WaitMessage = parfor_wait(size(labels,1), 'Waitbar', true, 'ReportInterval', 100);
            parfor i=1:size(labels,1) %length(obj.synSet(3).synLinks)
                tic
                % does the link that this synapse is on also have a synapse
                % in synSet2?
                [lia, locb] = ismember(obj.synSet(3).synLinks, obj.synSet(2).synLinks);
                % if yes, copy the ROI label of that synapse
                if lia
                    labels(i) = obj.synSet(2).synROIs(locb(1));
                else
                    % if still no luck, now you need to find this
                    % synapse's neighbors on synSet2. use synSet2
                    % because they will be closer. Don't have to search
                    % as far on the skeleton
                    neighbors = obj.synapseNeighbors(2, ...
                        obj.synSet(3).synLocs(i,:), ...
                        obj.synSet(3).synLinks(i), -1);
                    labels(i) = mode(obj.synSet(2).synROIs(neighbors(:,1)));
                end
                WaitMessage.Send;
                t=toc;
                fprintf('randsyn %d took %f s\n', i, t);
            end
            obj.synSet(3).synROIskey = obj.synSet(2).synROIskey;
            obj.synSet(3).synROIs = labels;
            WaitMessage.Destroy;
        end
        
        function obj = assignLinksToVoronoiSegments(obj, voronoiMask)
            % using a voronoiMask from a skeleton object (see the class
            % skeleton in the repository activityMap), assign every link and node in
            % the APL neurite skeleton to one of the Voronoi cells from the
            % 3d skeleton of the whole mushroom body
            numLinks = size(obj.links,1);
            numNodes = numLinks;
            obj.linkVoronoiSegments = zeros(numLinks,1);
            obj.nodeVoronoiSegments = zeros(numLinks,1);
            
            allCoords = zeros(length(obj.nodes),3);
            
            for i=1:length(obj.nodes)
                allCoords(i,:) = obj.nodes(i).coords;
            end
            
            rangeX = range(allCoords(:,1));
            rangeY = range(allCoords(:,2));
            rangeZ = range(allCoords(:,3));
            
            longestRange = max([rangeX, rangeY, rangeZ]);
            
            % squeeze everything so it fits in 256 x 256 pixels
            % this refers to code in connectomeToSkel.m
            
            divFactor = (longestRange+1)/256;
            
            minCoords = min(allCoords,[],1);
            
            for i=1:numLinks
                thisNode = obj.links(i,1);
                if thisNode % some of the link rows are all zeros, so skip them
                    if (obj.nodes(thisNode).coords(1) > 10000) & (obj.nodes(thisNode).coords(1) < 28000)
                        c = floor((obj.nodes(obj.links(i,1)).coords - minCoords)/divFactor)+1;
                        if c(3)<=size(voronoiMask,3)
                            voronoiSegment = find(voronoiMask(c(1),c(2),c(3),:));
                            if ~isempty(voronoiSegment)
                                obj.linkVoronoiSegments(i) = voronoiSegment;
                            end
                        end
                    end
                end
            end
            for i=1:numNodes
                if (obj.nodes(i).coords(1) > 10000) & (obj.nodes(i).coords(1) < 28000)
                    c = floor((obj.nodes(i).coords - minCoords)/divFactor) + 1;
                    if c(3)<=size(voronoiMask,3)
                        voronoiSegment = find(voronoiMask(c(1),c(2),c(3),:));
                        if ~isempty(voronoiSegment)
                            obj.nodeVoronoiSegments(i) = voronoiSegment;
                        end
                    end
                end
            end
            
        end
        
        function [meanRadius, pct5, pct95, maxRadius] = linkRadiusPerVoronoiSegment(obj)
            meanRadius = zeros(max(obj.linkVoronoiSegments),1);
            pct5 = meanRadius;
            pct95 = meanRadius;
%             stdRadius = meanRadius;
%             maxRadius = meanRadius;
%             minRadius = meanRadius;
            for i=1:length(meanRadius)
                theseLinks = (obj.linkVoronoiSegments==i);
                nodes1 = obj.links(theseLinks,1);
                nodes2 = obj.links(theseLinks,2);
                points1 = cell2mat({obj.nodes(nodes1).coords}');
                points2 = cell2mat({obj.nodes(nodes2).coords}');
                linkVectors = (points2 - points1);
                distances = sqrt(sum(linkVectors.^2,2));
                % radius of each segment is the mean of the radius of each
                % node at the ends
                radii = mean([obj.skel.radius(nodes1), obj.skel.radius(nodes2)], 2);
                meanRadius(i) = sum(radii.*distances)/sum(distances); %weighted mean
                pct5(i) = wprctile(radii,5,distances);
                pct95(i) = wprctile(radii,95,distances);
                maxRadius(i) = max(radii);
            end
        end
        
        function [normBranchPoints, normEndPoints] = normBranchPointsPerVoronoiSegment(obj)
            % number of branch points normalized to total length of
            % skeleton (in microns)
            normBranchPoints = zeros(max(obj.nodeVoronoiSegments),1);
            normEndPoints = normBranchPoints;
            linksPerNode = obj.linksPerNode();
            for i=1:length(normBranchPoints)
                theseLinks = (obj.linkVoronoiSegments==i); % OK to take the [0 0 0] links - they will have distance 0 and not add to the sum
                points1 = cell2mat({obj.nodes(obj.links(theseLinks,1)).coords}');
                points2 = cell2mat({obj.nodes(obj.links(theseLinks,2)).coords}');
                distances = sqrt(sum((points2 - points1).^2,2));
                totalLength = sum(distances);
                theseNodes = (obj.nodeVoronoiSegments==i);
                normBranchPoints(i) = sum((linksPerNode(theseNodes)-2).*(linksPerNode(theseNodes)>2)) / (totalLength*0.008);
                normEndPoints(i) = nnz(linksPerNode(theseNodes)==1) / (totalLength*0.008);
            end
        end
        
        function result = daughterLinks(obj)
            linksPerNode = obj.linksPerNode();
            branchPoints = find(linksPerNode>2);
            result = cell(length(branchPoints),1);
            for i=1:length(result)
                if mod(i,100)==0
                    i
                end
                neighbors = obj.nodes(branchPoints(i)).links;
                linkIndices = obj.getLinkIndices([neighbors', repmat(branchPoints(i),length(neighbors),1)]);
                radii = obj.links(linkIndices,3);
                result{i} = radii;
            end
        end
        
        function [resultAllreordered, totalPerm] = reorderKCs(obj,result)
            % operates on a numKCs x numKCs matrix called result
            % returns:
            % resultAllreordered: every square matrix is grouped by KC
            % subtype and each subtype is further grouped by hierarchical
            % clustering
            % totalPerm
            
            allNeurons = readtable('traced-neurons.csv');
            KCsubtypes = {
                'KCg';
                'KCg-d';
                'KCab-a';
                'KCab-c';
                'KCab-p';
                'KCab-s';
                'KCa''b''';
                };
            numTypes = length(KCsubtypes);
            KCsubtypeIndices = cell(numTypes,1);
            for i=1:numTypes
                %     KCsubtypes{i}
                [~,KCsubtypeIndices{i},~] = intersect(obj.synSet(1).uniquePartners, ...
                    allNeurons{strcmp(allNeurons{:,'instance'},KCsubtypes{i}),'bodyId'});
                %      length(KCsubtypeIndices{i})
                
            end
            
            % suppress automatic figure generation by dendrogram
            set(0,'DefaultFigureVisible','off');
            
            totalPerm = [];
            for i=1:numTypes
                
                resultThisType = result(KCsubtypeIndices{i},KCsubtypeIndices{i});
                %     size(resultThisType)
                z = linkage(resultThisType);
                h = figure;
                [~,~,perm] = dendrogram(z,length(KCsubtypeIndices{i}));
                %resultThisTypeReordered = resultThisType(perm,perm);
                totalPerm = [totalPerm; KCsubtypeIndices{i}(perm)];
                %figure,imagesc(resultThisTypeReordered);
                %     size(resultThisTypeReordered)
            end
            resultAllreordered = result(totalPerm,totalPerm);
            set(0,'DefaultFigureVisible','on');
            
        end
        function [resultAllreordered, totalPerm] = reorderKCsv11(obj,result)
            % operates on a numKCs x numKCs matrix called result
            % returns:
            % resultAllreordered: every square matrix is grouped by KC
            % subtype and each subtype is further grouped by hierarchical
            % clustering
            % totalPerm
            
            allNeurons = readtable('traced-neurons-v1.1.csv');
            KCsubtypes = {
                'KCg-m_R';
                'KCg-t_R';
                'KCg-d_R';
                'KCab-p_R';
                'KCab-s_R';
                'KCab-m_R';
                'KCab-c_R';
                'KCa''b''-ap1_R';
                'KCa''b''-ap2_R';
                'KCa''b''-m_R';
                };
            numTypes = length(KCsubtypes);
            KCsubtypeIndices = cell(numTypes,1);
            for i=1:numTypes
                %     KCsubtypes{i}
                [~,KCsubtypeIndices{i},~] = intersect(obj.synSet(1).uniquePartners, ...
                    allNeurons{strcmp(allNeurons{:,'instance'},KCsubtypes{i}),'bodyId'});
                %      length(KCsubtypeIndices{i})
                
            end
            
            % suppress automatic figure generation by dendrogram
            set(0,'DefaultFigureVisible','off');
            
            totalPerm = [];
            for i=1:numTypes
                
                resultThisType = result(KCsubtypeIndices{i},KCsubtypeIndices{i});
                     size(resultThisType)
                i
                z = linkage(resultThisType);
                h = figure;
                [~,~,perm] = dendrogram(z,length(KCsubtypeIndices{i}));
                %resultThisTypeReordered = resultThisType(perm,perm);
                totalPerm = [totalPerm; KCsubtypeIndices{i}(perm)];
                %figure,imagesc(resultThisTypeReordered);
                %     size(resultThisTypeReordered)
            end
            resultAllreordered = result(totalPerm,totalPerm);
            set(0,'DefaultFigureVisible','on');
            
        end        
        
    end
    
    methods (Static)
        function [invididuals, subtypes] = selfVsOthers(data, option, subtypeNumbers)
            % subtype numbers should be:
            %KCg	607
            %KCg-d	92
            %KCab-a	124
            %KCab-c	250
            %KCab-p	60
            %KCab-s	452
            %KCa'b'	336
            % give this as input: [607; 92; 124; 250; 60; 452; 336]
            
            if size(data,1)~=size(data,2)
                error('input data must be square matrix');
            end
            if sum(subtypeNumbers)~=size(data,1)
                error('sum of all subtypes must equal total number of neurons');
            end
            invididuals = zeros(size(data,1),1);
            subtypes = zeros(length(subtypeNumbers),1);
            runningCount = 0;
            for j=1:length(subtypeNumbers)
                begini = runningCount + 1;
                endi = runningCount + subtypeNumbers(j);
                runningCount = endi;
                thisSubtype = data(begini:endi,begini:endi);
                meanThisSubtype = mean(thisSubtype(:));
                switch option
                    case 'rows'
                        for i=begini:endi
                            self = data(i,i);
                            others = mean(data(i,setdiff(begini:endi,i)));
                            invididuals(i) = self/others;
                        end
                        otherSubtype = data(begini:endi,setdiff(1:end,begini:endi));
                        subtypes(j) = meanThisSubtype/mean(otherSubtype,'all');
                    case 'columns'
                        for i=begini:endi
                            self = data(i,i);
                            others = mean(data(setdiff(begini:endi,i),i));
                            invididuals(i) = self/others;
                        end
                        otherSubtype = data(setdiff(1:end,begini:endi),begini:endi);
                        subtypes(j) = meanThisSubtype/mean(otherSubtype,'all');
                    otherwise
                        error('invalid option');
                end
            end
        end
    end
    
    
end