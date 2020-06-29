function result = convolveExpDecay(input, skel2d, decayFactor)
% convolveExpDecay(input, skel2d, decayFactor)
% input: red dye stimulation (vstimFit, etc)
% skel2d: coordinates of the backbone skeleton in 2 dimensions in microns
% decayFactor: space constant in microns

result = zeros(length(input),1);
for i=1:length(input)
    result = result + input(i)*exp(-pdist2(skel2d,skel2d(i,:),'cityblock')/decayFactor);
end
end