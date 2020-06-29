% this script operates on a numKCs x numKCs matrix called result

allNeurons = readtable('../hemibrain-synapse-spacing/data/traced-neurons.csv');
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
    [~,KCsubtypeIndices{i},~] = intersect(apl.synSet(1).uniquePartners, ...
        allNeurons{strcmp(allNeurons{:,'instance'},KCsubtypes{i}),'bodyId'});
%      length(KCsubtypeIndices{i})

end
meanResult = zeros(numTypes,2,size(result,3));
resultAllreordered = cell(size(result,3),1);
set(0,'DefaultFigureVisible','off');

for j=1:size(result,3)
    totalPerm = [];
    for i=1:numTypes
        
        resultThisType = result(KCsubtypeIndices{i},KCsubtypeIndices{i},j);
        meanResult(i,1,j) = mean(resultThisType(eye(size(resultThisType))==1));
        meanResult(i,2,j) = mean(resultThisType(eye(size(resultThisType))==0));
        %     size(resultThisType)
        z = linkage(resultThisType);
        h = figure;
        [h,t,perm] = dendrogram(z,length(KCsubtypeIndices{i}));
        resultThisTypeReordered = resultThisType(perm,perm);
        totalPerm = [totalPerm; KCsubtypeIndices{i}(perm)];
        %figure,imagesc(resultThisTypeReordered);
        %     size(resultThisTypeReordered)
    end
    resultAllreordered{j} = result(totalPerm,totalPerm,j);
end
set(0,'DefaultFigureVisible','on');


