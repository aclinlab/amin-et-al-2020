allCoords = zeros(length(apl.nodes),3);

for i=1:length(apl.nodes)
    allCoords(i,:) = apl.nodes(i).coords;
end

rangeX = range(allCoords(:,1));
rangeY = range(allCoords(:,2));
rangeZ = range(allCoords(:,3));

longestRange = max([rangeX, rangeY, rangeZ]);

% squeeze everything so it fits in 256 x 256 pixels

divFactor = (longestRange+1)/256;

minCoords = min(allCoords,[],1);

allCoords = allCoords - repmat(minCoords,size(allCoords,1),1);

allCoords = floor(allCoords/divFactor) + 1;

maskSize = ceil([rangeX, rangeY, rangeZ]/divFactor);

mask = false(maskSize);

for i=1:length(apl.nodes)
    xyz = floor((apl.nodes(i).coords - minCoords)/divFactor) + 1;
    xyz = xyz + repmat([-1; 0; 1],3,1);
    xyz(xyz<1) = 1;
    xyz(xyz>256) = 256;
    mask(xyz(:,1),xyz(:,2),xyz(:,3)) = true;
end
for i=1:length(apl.synSet(2).synLinks)
    xyz = floor((apl.synSet(2).synLocs(i,:) - minCoords)/divFactor) + 1;
    xyz = xyz + repmat([-1; 0; 1],3,1);
    xyz(xyz<1) = 1;
    xyz(xyz>256) = 256;
    mask(xyz(:,1),xyz(:,2),xyz(:,3)) = true;
end

pixelCal = repmat(0.008*divFactor,1,3); % each pixel before was 8 nm but now we are binning together divFactor pixels in each pixel

% ignore the bottom 96 pixels - there is no mushroom body there
aplSkel = skeleton(mask(:,:,1:160),pixelCal);
% this is where you have to draw the skeleton manually

nodeLinks = {aplSkel.nodes.links};
numLinksPerNode = zeros(length(nodeLinks),1);
for j=1:length(nodeLinks)
    numLinksPerNode(j) = length(nodeLinks{j});
end

[maxLinks,junctionIndex] = max(numLinksPerNode);
if maxLinks~=3
    error('the node with the most links does not have 3 links!');
end

% register the connectome APL onto our standard average mushroom body 3D
% skeleton
aplEndNodes = aplSkel.findEndNodes();

[aplDistances,~] = aplSkel.getDistances(junctionIndex);

aplLengths = aplDistances(aplEndNodes);

normalizedLengths = [97.61748775; 181.6842957; 73.20134649]; % these are the lengths from our standard average skeleton
stretchFactors = aplLengths./normalizedLengths;
aplSkel = aplSkel.labelBranches;
% Re-draw the  spaced nodes, emanating from the junction
% Set the starting point to be the junction
aplSkel.userStartingPoint = aplSkel.nodes(junctionIndex).realCoords;
aplSkel = aplSkel.createSpacedNodes(10*stretchFactors);
aplSkel = aplSkel.createVoronoiMask();

