% gmmem.m

%
% Gaussian Mixture Model using Expectation Maximisation
%

function [p, m, C] = gmmem(data, p0, m0, C0)

% Constants

[nDimensions, nData] = size(data);
nClasses             = length(p0);

% minDetCov = max(1E-20, 10^(-200/nDimensions)); % minimum determinant of covariance matrix
minDetCov = 0; 

if (nClasses <= 1)
  % Compute gaussian
  p(1) = 1;
  m{1} = mean(data')';
  C{1} = cov(data');
  return;
end;

%
% Compute responsibilities p( class | data )
%

% for i = 1:nClasses
%   CInvRoot = chol(inv(C0{i}));
%   CDetRoot = sqrt(det(C0{i}));
%   cData    = data - m0{i} * ones(1, nData);
%   ei = CInvRoot * cData;
%   ei = 0.5 * sum(ei .^2);
%   respi = p0(i) * CDetRoot * exp( - ei );
%   resp(i, :) = respi;
% end;

% TEST 
for i = 1:nClasses
  CInvRoot = chol(inv(C0{i}));
  LogDetC  = log(det(C0{i}));
  cData    = data - m0{i} * ones(1, nData);
  ei = CInvRoot * cData;
  ei = ei .* ei;
  ei = sum(ei);
  % ei = ei ./ (nDimensions/3);
  ei = ei + LogDetC * ones(1, nData);
  ei = ei - 2 * (log(p0(i)) * ones(1, nData));
  resp(i, :) = exp(-ei);
end;

% test for valid points (responsibility > 0)
validPoint = sum(resp, 1) > 0;
data = data(:, validPoint);
resp = resp(:, validPoint);
[nDimensions, nData] = size(data);

resp = resp ./ (ones(nClasses, 1) * sum(resp, 1));

%
% Compute mean, covariance and priors
%

sumresp = sum(resp, 2);
p = (sumresp / nData)';

validClass = ones(nClasses, 1);
for i = 1:nClasses
  if (p(i) <= 0)
    validClass(i) = 0;
    continue;
  end;
  
  respi = ones(nDimensions, 1) * resp(i, :);
  m{i}  = sum(respi .* data, 2) / sumresp(i);
  cData = data - m{i} * ones(1, nData);
  C{i}  = (respi .* cData) * cData' / sumresp(i);
  
  % test for valid covariance
  if (det(C{i}) < minDetCov)
    validClass(i) = 0;
  else
    [CInvRoot, err] = chol(inv(C{i}));
    if (err)
      validClass(i) = 0;
    end;
  end;

end;

validClassIndex = find(validClass);
p = p(:, validClassIndex);
m = m(:, validClassIndex);
C = C(:, validClassIndex);








