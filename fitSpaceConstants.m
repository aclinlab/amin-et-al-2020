% these data are also stored in Fig 8.xlsx associated with the paper
load('redDyeVoronoiFits.mat'); %vstimFit, hstimFit, cstimFit
load('voronoiIndices.mat'); %ctohindices, ctovindices
load('skel2d.mat'); %skel2d
redDyeFits = [hstimFit, vstimFit, cstimFit];
outputNames = {'hstim', 'vstim', 'cstim'};
% these are values of k in units of sqrt(pixels)
% to convert to sqrt(m), multiply by sqrt(8e-9 m/pixel) (8 nm = 1 pixel)
electroDistConstants = [625, 1250, 1875, 2381, 11180]; % k = .0559, 0.1118, 0.1677, 0.213, 1.0 m
% the value k=0.213 m comes from Cuntz et al 2013 where Rm = 0.8166 Ohm*m^2
% and Ra = 9 Ohm*m. k=1.0 m comes from Gouwens et al 2009, Rm = 2.04 Ohm*m^2
% and Ra = 1.02 Ohm*m.
realDistConstants = electroDistConstants*5; %this gives real distance space constant in pixels
% to convert to m, multiply by 8e-9 m/pixel
% 25, 50, 75, 95.2, 223.6, 447.2 um
micronConstants = realDistConstants*8e-9/1e-6; % space constants in microns
EDresults = nan(length(ctovindices),2*length(electroDistConstants), size(redDyeFits,2));
realDistResults = EDresults;
backboneResult = EDresults;
% loop over vstim, hstim, cstim
for i=1:size(redDyeFits,2)
    for j=1:length(electroDistConstants)
        % electrotonic distance
        [result,~]=apl.stimVoronoi(electroDistConstants(j),randDistED,3,redDyeFits(:,i));
        EDresults(1:length(ctohindices),(2*j)-1,i) = result(ctohindices)/max(result);
        EDresults(1:length(ctovindices),(2*j),i) = result(ctovindices)/max(result);
        % real distance
        [result,~]=apl.stimVoronoi(realDistConstants(j),randDistReal,3,redDyeFits(:,i));
        realDistResults(1:length(ctohindices),(2*j)-1,i) = result(ctohindices)/max(result);
        realDistResults(1:length(ctovindices),(2*j),i) = result(ctovindices)/max(result);
        % distance on backbone skeleton
        result = convolveExpDecay(redDyeFits(:,i),skel2d,micronConstants(j));
        backboneResult(1:length(ctohindices),(2*j)-1,i) = result(ctohindices)/max(result);
        backboneResult(1:length(ctovindices),(2*j),i) = result(ctovindices)/max(result);
        
    end
    csvwrite(strcat(outputNames{i},'-','ED',datestr(now,'yyyymmddTHHMM'),'.csv'), EDresults(:,:,i));
    csvwrite(strcat(outputNames{i},'-','real',datestr(now,'yyyymmddTHHMM'),'.csv'), realDistResults(:,:,i));
    csvwrite(strcat(outputNames{i},'-','backbone',datestr(now,'yyyymmddTHHMM'),'.csv'), backboneResult(:,:,i));
end
