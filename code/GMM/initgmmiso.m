% initgmmiso.m

% Initialise mixture components evenly distributed along 
% 1st principal direction

function [p, m, sigma] = initgmmiso(data, nClasses, options)

[nDimensions, nData] = size(data);

if (nClasses == 1)
  % Compute isotropic gaussian
  p(1) = 1;
  m{1} = mean(data')';
  cData = data - m{1} * ones(1, nData);
  sumsq = sum(sum(cData .^ 2));
  sigma(1) = sqrt(sumsq / (nData * nDimensions));
  return;
end;

CData = cov(data');
[V, D] = eig(CData);

v1 = V(:, nDimensions);
d1 = D(nDimensions, nDimensions);

step         = nData / nClasses;
stepIndexTmp = step / 2;
for i = 1:nClasses
  stepIndex(i) = ceil(stepIndexTmp);
  stepIndexTmp = stepIndexTmp + step;
end;

[sortedData, sortIndex] = sort(v1' * data);

for i = 1:nClasses
  p(i)     = 1 / nClasses;
  if (options.init == 0)
     m{i} = data(:, sortIndex(stepIndex(i))); % this is the middle data value of the nth classe
  else
     m{i} = rand(nDimensions,1);
  end
  sigma(i) = sqrt(d1) / nClasses; 
end;
