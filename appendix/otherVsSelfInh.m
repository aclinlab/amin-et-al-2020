% This creates panel C-E of Appendix Figure 1 in Amin et al., 2020

numRand = 100; % number of neurons
sims = 0.1:0.1:1; % number is similarity values (sigma in the figure) to test
thetas = 0:0.05:1; % number of thresholds to test
alpha = 0.5; % gain on inhibition
numTrials = 100; % number of simulations to carry out
[noinhRand, selfinhRand, allinhRand, noinhCL, selfinhCL, allinhCL, noinhSP, selfinhSP, allinhSP] = ...
    deal(zeros(length(thetas),length(sims), numTrials));
% selfinhRand = zeros(length(thetas),length(sims));
% allinhRand = zeros(length(thetas),length(sims));
% noinhCL = zeros(length(thetas),length(sims));
% selfinhCL = zeros(length(thetas),length(thetas));
% allinhCL = zeros(length(thetas),length(thetas));
for k=1:numTrials
for i=1:length(thetas)
    for j=1:length(sims)
        % create the random neuron activities. x is x(1) in the figure
        x = rand(1, numRand);
        y = x+sims(j)*randn(1, numRand); % y is x(2) in the figure
        exc = [x; y];
        
        % this is the activity of the neurons with no inhibition
        noinhAct = exc-thetas(i);
        noinhAct(noinhAct<0) = 0;
        noinhRand(i,j,k) = pdist(noinhAct,'cosine'); % cosine distances
        noinhCL(i,j,k) = nnz(noinhAct)/numel(noinhAct); % coding level (fraction of active neurons)
        
        % activity of neurons with self-inhibition
        selfinhAct = (1-alpha)*exc - thetas(i);
        selfinhAct(selfinhAct<0) = 0;
        selfinhRand(i,j,k) = pdist(selfinhAct, 'cosine');
        selfinhCL(i,j,k) = nnz(selfinhAct)/numel(selfinhAct);
        
        % activity of neurons with other-inhibition
        allinhAct = exc - repmat(mean(exc,2),1,size(exc,2))*alpha - thetas(i);
        allinhAct(allinhAct<0) = 0;
        allinhRand(i,j,k) = pdist(allinhAct, 'cosine');
        allinhCL(i,j,k) = nnz(allinhAct)/numel(allinhAct);
    end
end
end

meanNoinhRand = nanmean(noinhRand,3);
meanSelfinhRand = nanmean(selfinhRand,3);
meanAllinhRand = nanmean(allinhRand,3);
meanNoinhCL = nanmean(noinhCL,3);
meanSelfinhCL = nanmean(selfinhCL,3);
meanAllinhCL = nanmean(allinhCL,3);
% panel C
figure, 
subplot(1,3,1), imagesc(meanNoinhRand, [0 1]), axis equal, axis tight;
subplot(1,3,2), imagesc(meanSelfinhRand, [0 1]), axis equal, axis tight;
subplot(1,3,3), imagesc(meanAllinhRand, [0 1]), axis equal, axis tight;
% panel D
figure, 
subplot(1,3,1), imagesc(meanNoinhCL, [0 1]), axis equal, axis tight;
subplot(1,3,2), imagesc(meanSelfinhCL, [0 1]), axis equal, axis tight;
subplot(1,3,3), imagesc(meanAllinhCL, [0 1]), axis equal, axis tight;

% panel E
figure,
hold on
for i=1:length(sims)
    plot(meanNoinhCL(:,i),meanNoinhRand(:,i),'o-b');
    plot(meanSelfinhCL(:,i),meanSelfinhRand(:,i),'o-r');
    plot(meanAllinhCL(:,i),meanAllinhRand(:,i),'o-g');
end

