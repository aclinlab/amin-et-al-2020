% This creates panel B of appendix figure 1 of Amin et al., 2020
num = 20; % different values of x_0 in the figure
js = 0.1:0.1:0.9; % different values of a in the figure
noinh = zeros(num, length(js));
selfinh = zeros(num, length(js));
allinh = zeros(num, length(js));
noinhActAvg = zeros(num, length(js));
selfinhActAvg = zeros(num, length(js));
allinhActAvg = zeros(num, length(js));

alpha = 0.5;
theta = 2;

for i=1:num
    for j=1:length(js)
        % excitation onto neurons x(1) and x(2) in the figure
        exc = [i, i*js(j); 
               i*js(j), i];
        noinhAct = exc-theta;
        noinhAct(noinhAct<0) = 0;
        selfinhAct = (1-alpha)*exc - theta;
        selfinhAct(selfinhAct<0) = 0;
        allinhAct = exc - mean([i, i-js(j)])*alpha - theta;
        allinhAct(allinhAct<0) = 0;
        noinh(i,j) = pdist(noinhAct, 'cosine');
        selfinh(i,j) = pdist(selfinhAct, 'cosine');
        allinh(i,j) = pdist(allinhAct, 'cosine');
    end
end

% panel B of the figure
figure,
subplot(1,3,1), imagesc(flipud(noinh)), axis equal, axis tight;
subplot(1,3,2), imagesc(flipud(selfinh)), axis equal, axis tight;
subplot(1,3,3), imagesc(flipud(allinh)), axis equal, axis tight;
