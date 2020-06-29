function [] = writeTiffParulaWithRange(data, filename, range)

data = double(data);
% crop the limits
data(data>range(2)) = range(2);
data(data<range(1)) = range(1);

% normalize the data to lie betwen 0 and 1
data = (data-range(1))/(range(2)-range(1));

imwrite(data*64, parula, ...
    filename, ...
    'tif', 'Compression', 'lzw');

end