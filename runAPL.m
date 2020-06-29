% NB do not just run this code - it will take a
% very long time and also several sections are interactive (need to copy
% and paste data from a spreadsheet)
% Rather provided as a reference to how to produce the data
% in Amin et al. Code running times from a 2013 MacBook Pro.

% load the skeleton (few minutes)
apl = APLskel();
% analyse skeleton for unconnected fragments
apl = apl.labelNodes(); 

% find where all the synapses are on the skeleton (~15-30 minutes)
apl = apl.getKCSynapses();
apl = apl.labelSynsByROI([1 2]);

% Numbers from the text:
totalLength = apl.getTotalLengthMainBranch(0);
totalEDLength = apl.getTotalLengthMainBranch(1);

% add 'random' locations on APL
apl = apl.addEvenRandomChans(10000, 3);

% Figure 8B-D:
apl.drawSparseSynapses(2000, 1, 'r'); % APL-KC synapses
apl.drawSparseSynapses(2000, 2, 'b'); % KC-APL synapses
apl.drawSparseSynapsesAPL(2000, 3); % 'random' locations on APL

% Process the random synapses
apl = apl.getSynSetGraph(3); % ~3 h for 10,000 synapses
save(strcat('apl',datestr(now,'yyyymmddTHHMMSS'),'.mat'),'apl');
randDistED = apl.getAllSynDistWithinSynSet(3,1e7,3); % ~30' for 10k synapses
save(strcat('aplRandDistED',datestr(now,'yyyymmddTHHMMSS'),'.mat'),'randDistED');
randDistReal = apl.getAllSynDistWithinSynSet(3,1e7,2);
save(strcat('aplRandDistReal',datestr(now,'yyyymmddTHHMMSS'),'.mat'),'randDistReal');
apl = apl.labelRandomSynsByROI(); % ~70' for 10k synapses
save(strcat('apl',datestr(now,'yyyymmddTHHMMSS'),'.mat'),'apl');

% map the connectome APL onto our standard mushroom body 3D skeleton
connectomeToSkel; % or load aplmanualSkel200303.mat --> aplSkel

% assign nodes and links of APL skeleton to Voronoi segments according to
% the manual skeleton from connectomeToSkel
apl = apl.assignLinksToVoronoiSegments(aplSkel.voronoiMask);
apl = apl.assignSynapsesToVoronoiSegments();

% Fig 8F
fitSpaceConstants;

% make a graph of each set of synapses (APLtoKC, KCtoAPL) with distances
% between all of them along the skeleton (takes ~7.5 h for both synSets)
apl = apl.getSynSetGraph(1); save(strcat('apl',datestr(now,'yyyymmddTHHMMSS'),'.mat'),'apl'); % 5 hours
apl = apl.getSynSetGraph(2); save(strcat('apl',datestr(now,'yyyymmddTHHMMSS'),'.mat'),'apl'); % 2.5 hours
% use Dijkstra's algorithm to get distances between every synapse in the
% APLtoKC set (~3 h). The file will be ~865 MB
EDMatrix = apl.getAllSynDistWithinSynSet(1,1e7,3); 
save(strcat('APLtoKCEDMatrix',datestr(now,'yyyymmddTHHMMSS'),'.mat'), 'EDMatrix');

% for each KCtoAPL synapse, find its nearest neighbors in the APLtoKC set.
% (~18 h). The file will be ~12 MB
KCtoAPLneighborsOnAPLtoKC = apl.getNeighborsOnOtherGraph([2 1]);
save(strcat('KCtoAPLneighborsOnAPLtoKC',datestr(now,'yyyymmddTHHMMSS'),'.mat'),...
    'KCtoAPLneighborsOnAPLtoKC');



% combine distMatric and KCtoAPLneighborsOnAPLtoKC to find distance from
% each KCtoAPL synapse to each APltoKC synapse (this is fast: 6-7 min)
KCtoAPLEDtoAPLtoKC = apl.getDistBetweenSynSets(EDMatrix,KCtoAPLneighborsOnAPLtoKC,1); 
% need to use -v.3 because this is a huge file: ~5.3 GB
save(strcat('KCtoAPLEDtoAPLtoKC',datestr(now,'yyyymmddTHHMMSS'),'.mat'), ...
        'KCtoAPLEDtoAPLtoKC','-v7.3');

apl = apl.mapOrigSynsToSynLocs;

% for this one, load KCtoAPLEDtoAPLtoKC into memory
runWeightedSumBetweenSynSets;

realDistMatrix = apl.getAllSynDistWithinSynSet(1,1e7,2); 
save(strcat('APLtoKCrealDistMatrix',datestr(now,'yyyymmddTHHMMSS'),'.mat'), 'realDistMatrix');
% real distance
KCtoAPLrealDistToAPLtoKC = apl.getDistBetweenSynSets(realDistMatrix,KCtoAPLneighborsOnAPLtoKC,0); 
% need to use -v.3 because this is a huge file: ~5.3 GB
save(strcat('KCtoAPLrealDistToAPLtoKC',datestr(now,'yyyymmddTHHMMSS'),'.mat'), ...
        'KCtoAPLrealDistToAPLtoKC','-v7.3');

% for this one, load KCtoAPLrealDistToAPLtoKC into memory
runWeightedSumBetweenSynSetsRealDist;

% Figure 8-fig supp 3:
[meanRadius, pct5, pct95, maxRadius] = apl.linkRadiusPerVoronoiSegment;