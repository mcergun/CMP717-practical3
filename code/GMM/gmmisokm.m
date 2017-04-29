% gmmisokm.m

%
% Gaussian Mixture Model with Isotropic Covariance Matrices 
% using K-means
%

function [p, m, sigma] = gmmisokm(data, p0, m0, sigma0)

% Constants

[nDimensions, nData] = size(data);
nClasses             = length(p0);

minSigma = max(1E-10, 10^(-100/nDimensions)); % minimum standard deviation

if (nClasses <= 1)
  % Compute isotropic gaussian
  p(1) = 1;
  m{1} = mean(data')';
  cData = data - m{1} * ones(1, nData);
  sumsq = sum(sum(cData .^ 2));
  sigma(1) = sqrt(sumsq / (nData * nDimensions));
  return;
end;

% Compute error e
% e = - log( p(class | data) ) + const

for i = 1:nClasses
  cData = data - m0{i} * ones(1, nData);
  ei = cData ./ sigma0(i);
  ei = 0.5 * sum(ei .^2);
  ei = ei + nDimensions * log( sigma0(i) );
  ei = ei - log(p0(i));
  e(i, :) = ei;
end;

% Assign class with mimimum error to each point
% and sort the classes

[minError, class]  = min(e);
[classSort, order] = sort(class);
classSortShiftR(2:nData + 1) = classSort(1:nData);
classSortShiftR(1) = 0;
classStart = find( classSortShiftR(1:nData) ~= classSort );

% Recompute class mean, variance and prior

nClasses = length(classStart);
classStart(nClasses + 1) = nData + 1;

validClass = ones(1, nClasses);
for i = 1:nClasses
  classStarti   = classStart(i);
  classStartip1 = classStart(i + 1);
  nClass        = classStartip1 - classStarti;
  
  if (nClass < nDimensions) 
    % covariance degenerate
    validClass(1, i) = 0;
    continue;
  end;
  
  orderi        = order(classStarti : classStartip1 - 1);
  datai         = data(:, orderi);

  m{i}      = mean( datai' )';  
  cDatai    = datai - m{i} * ones(1, nClass);
  sumsq     = sum(sum(cDatai .* cDatai));
  sigma(i)  = sqrt(sumsq / (nClass * nDimensions));
  p(i)  = (classStartip1 - classStarti) / nData;
  
  if (sigma(i) < minSigma)
    validClass(1, i) = 0;
  end;

end;

validClassIndex = find(validClass);
p      = p   (:, validClassIndex);
m      = m   (:, validClassIndex);
sigma = sigma(:, validClassIndex);

