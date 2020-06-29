electroDistConstants = [625, 1250, 1875, 2381, 5590, 11180, Inf];
realDistConstants = electroDistConstants*5;

% This section contains examples of commands used to generae Fig 8I-L (but
% not every condition)
% Fig 8J, K
numKCs = 1929;

resultAll = zeros(numKCs, numKCs, length(electroDistConstants));
% This takes about 90'
for i=1:length(electroDistConstants)
    resultAll(:,:,i) = apl.weightedSumBetweenSynSets(electroDistConstants(i), KCtoAPLEDtoAPLtoKC, [2 1],...
        {{'CA(R)','PED(R)','gL(R)','aL(R)', 'bL(R)', 'a''L(R)', 'b''L(R)'}});
end
save(strcat('results_allKCs_',datestr(now,'yyyymmddTHHMMSS'),'.mat'),'resultAll');

% the following script is used to sort KCs by subtype:
[resultsAllreordered, totalPerm1250] = apl.reorderKCs(resultAll(:,:,2)); % use 1250
% save the order of KCs from 50 um, all KCs (Fig 8G)

% Fig 8G
writeTiffParulaWithRange(resultAll(totalPerm1250,totalPerm1250,2), 'allKCs.tiff',[0 500]);
% Fig 8H
writeTiffParulaWithRange(resultAll(totalPerm1250,totalPerm1250,7), 'allKCsNoDecay.tiff',[0 8000]);

subtypeNumbers = [607 %g
    92 %g-d
    124 %ab-a
    250 %ab-c
    60 %ab-p
    452 %ab-s
    336]; %a'b'

% Fig 8J - self vs all KCs
% Fig 8K - self vs only KCs of the same subtype
selfVsAllKCs = zeros(length(totalPerm1250), length(electroDistConstants)+1);
selfVsSameTypeKCs = zeros(length(totalPerm1250), length(electroDistConstants)+1);
for i=1:length(electroDistConstants)
    [selfVsAllKCs(:,i), ~] = apl.selfVsOthers(resultAll(totalPerm1250,totalPerm1250,i),'columns',length(totalPerm1250));
    [selfVsSameTypeKCs(:,i), ~] = apl.selfVsOthers(resultAll(totalPerm1250,totalPerm1250,i),'columns',subtypeNumbers);
end

aplScramble = apl.scrambleUniquePartners([1 2]);
resultScrambled = aplScramble.weightedSumBetweenSynSets(625, KCtoAPLEDtoAPLtoKC, [2 1],...
    {{'CA(R)','PED(R)','gL(R)','aL(R)', 'bL(R)', 'a''L(R)', 'b''L(R)'}});
save('result_KCsScrambled_lambda625.mat','resultScrambled');
[selfVsAllKCs(:,length(electroDistConstants)+1), ~] = ...
    apl.selfVsOthers(resultScrambled(totalPerm1250,totalPerm1250),'columns',length(totalPerm1250));
[selfVsSameTypeKCs(:,length(electroDistConstants)+1), ~] = ...
    apl.selfVsOthers(resultScrambled(totalPerm1250,totalPerm1250),'columns',subtypeNumbers);

% Fig 8J - self vs all KCs in the order of electroDistConstants with
% scrambled-KC-indices last
csvwrite(strcat('selfVsAllKCs_',datestr(now,'yyyymmddTHHMM'),'.csv'), selfVsAllKCs);
% Fig 8K
csvwrite(strcat('selfVsSameTypeKCs_',datestr(now,'yyyymmddTHHMM'),'.csv'), selfVsSameTypeKCs);

resultCalyx = zeros(numKCs, numKCs, length(electroDistConstants));
% each loop will take ~15-20'
for i=1:length(electroDistConstants)
    resultCalyx(:,:,i) = apl.weightedSumBetweenSynSets(electroDistConstants(i), KCtoAPLEDtoAPLtoKC, [2 1],...
        {'CA(R)'});
end
save(strcat('results_calyx_',datestr(now,'yyyymmddTHHMMSS'),'.mat'),'resultCalyx');
resultCalyxScrambled = aplScramble.weightedSumBetweenSynSets(625, KCtoAPLEDtoAPLtoKC, [2 1],...
    {'CA(R)'});
save('result_KCsScrambled_lambda625_calyx.mat','resultCalyxScrambled');
    
% Fig 8I
writeTiffParulaWithRange(resultCalyx(totalPerm1250,totalPerm1250,2), 'calyx1250.tiff',[0 200]);

% Fig 8L
selfVsSameTypeKCsCalyx = zeros(length(totalPerm1250), length(electroDistConstants)+1);
for i=1:length(electroDistConstants)
    [selfVsSameTypeKCsCalyx(:,i), ~] = apl.selfVsOthers(resultCalyx(totalPerm1250,totalPerm1250,i),'columns',subtypeNumbers);
end
[selfVsSameTypeKCsCalyx(:,length(electroDistConstants)+1), ~] = ...
    apl.selfVsOthers(resultCalyxScrambled(totalPerm1250,totalPerm1250),'columns',subtypeNumbers);
csvwrite(strcat('selfVsSameTypeKCsCalyx_',datestr(now,'yyyymmddTHHMM'),'.csv'), selfVsSameTypeKCsCalyx);



