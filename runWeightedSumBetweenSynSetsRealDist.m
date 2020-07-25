electroDistConstants = [624, 1249, 1873, 2381, 5590, 11180, Inf];
realDistConstants = electroDistConstants*5.0062;

% This section contains examples of commands used to generae Fig 8I-L (but
% not every condition)
% Fig 8J, K
numKCs = length(apl.synSet(1).uniquePartners);

resultAll = zeros(numKCs, numKCs, length(realDistConstants));
% eThis takes about 90'
for i=1:length(realDistConstants)
    resultAll(:,:,i) = apl.weightedSumBetweenSynSets(realDistConstants(i), KCtoAPLrealDistToAPLtoKC, [2 1],...
        {apl.synSet(1).synROIskey}); %{'CA(R)','PED(R)','gL(R)','aL(R)', 'bL(R)', 'a''L(R)', 'b''L(R)'}
end
save(strcat('results_allKCs_realDist_',datestr(now,'yyyymmddTHHMMSS'),'.mat'),'resultAll');

% the following script is used to sort KCs by subtype:
[resultsAllreordered, totalPerm6250] = apl.reorderKCsv11(resultAll(:,:,2)); % use 6250
% save the order of KCs from 50 um, all KCs (Fig 8G)

% Fig 8G
writeTiffParulaWithRange(resultAll(totalPerm6250,totalPerm6250,2), strcat('allKCs_realDist-',datestr(now,'yyyymmddTHHMM'),'.tiff'),[0 500]);
% Fig 8H
writeTiffParulaWithRange(resultAll(totalPerm6250,totalPerm6250,6), strcat('allKCsNoDecay_realDist-',datestr(now,'yyyymmddTHHMM'),'.tiff'),[0 8000]);

subtypeNumbers = [590 %g-m
    8 %g-t
    99 %g-d
    60 %ab-p
    223 %ab-s
    354 %ab-m
    252 %ab-c
    91 %a'b'-ap1
    127 %a'b'-ap2
    119]; %a'b'-m

% Fig 8J - self vs all KCs
% Fig 8K - self vs only KCs of the same subtype
selfVsAllKCs = zeros(length(totalPerm6250), length(realDistConstants)+1);
selfVsSameTypeKCs = zeros(length(totalPerm6250), length(realDistConstants)+1);
for i=1:length(realDistConstants)
    [selfVsAllKCs(:,i), ~] = apl.selfVsOthers(resultAll(totalPerm6250,totalPerm6250,i),'columns',length(totalPerm6250));
    [selfVsSameTypeKCs(:,i), ~] = apl.selfVsOthers(resultAll(totalPerm6250,totalPerm6250,i),'columns',subtypeNumbers);
end

aplScramble = apl.scrambleUniquePartners([1 2]);
resultScrambled = aplScramble.weightedSumBetweenSynSets(3125, KCtoAPLrealDistToAPLtoKC, [2 1],...
    {apl.synSet(1).synROIskey});
save('result_KCsScrambled_realDist_lambda3125.mat','resultScrambled');
[selfVsAllKCs(:,length(realDistConstants)+1), ~] = ...
    apl.selfVsOthers(resultScrambled(totalPerm6250,totalPerm6250),'columns',length(totalPerm6250));
[selfVsSameTypeKCs(:,length(realDistConstants)+1), ~] = ...
    apl.selfVsOthers(resultScrambled(totalPerm6250,totalPerm6250),'columns',subtypeNumbers);

% Fig 8J - self vs all KCs in the order of realDistConstants with
% scrambled-KC-indices last
csvwrite(strcat('selfVsAllKCs_realDist_',datestr(now,'yyyymmddTHHMM'),'.csv'), selfVsAllKCs);
% Fig 8K
csvwrite(strcat('selfVsSameTypeKCs_realDist_',datestr(now,'yyyymmddTHHMM'),'.csv'), selfVsSameTypeKCs);

resultCalyx = zeros(numKCs, numKCs, length(realDistConstants));
% each loop will take ~15-20'
for i=1:length(realDistConstants)
    resultCalyx(:,:,i) = apl.weightedSumBetweenSynSets(realDistConstants(i), KCtoAPLrealDistToAPLtoKC, [2 1],...
        {{'CA(R)','','SLP(R)','SCL(R)','ICL(R)','PLP(R)'}},1);
end
save(strcat('results_calyx_',datestr(now,'yyyymmddTHHMMSS'),'.mat'),'resultCalyx');
resultCalyxScrambled = aplScramble.weightedSumBetweenSynSets(3125, KCtoAPLrealDistToAPLtoKC, [2 1],...
    {{'CA(R)','','SLP(R)','SCL(R)','ICL(R)','PLP(R)'}},1);
save('result_KCsScrambled_realDist_lambda3125_calyx','resultCalyxScrambled');
    
% Fig 8I
writeTiffParulaWithRange(resultCalyx(totalPerm6250,totalPerm6250,2), strcat('calyx6250_realDist-',datestr(now,'yyyymmddTHHMM'),'.tiff'),[0 200]);

% Fig 8L
selfVsSameTypeKCsCalyx = zeros(length(totalPerm6250), length(realDistConstants)+1);
for i=1:length(realDistConstants)
    [selfVsSameTypeKCsCalyx(:,i), ~] = apl.selfVsOthers(resultCalyx(totalPerm6250,totalPerm6250,i),'columns',subtypeNumbers);
end
[selfVsSameTypeKCsCalyx(:,length(realDistConstants)+1), ~] = ...
    apl.selfVsOthers(resultCalyxScrambled(totalPerm6250,totalPerm6250),'columns',subtypeNumbers);
csvwrite(strcat('selfVsSameTypeKCsCalyx_realDist_',datestr(now,'yyyymmddTHHMM'),'.csv'), selfVsSameTypeKCsCalyx);



