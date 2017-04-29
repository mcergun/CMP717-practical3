% DrawEllipse3D.m

%
% Draws the 3D ellipse
%
% (x - m) * C-1 * (x - m)' = 0 
%

function DrawEllipse3D( m, C, colour );

if (nargin < 3)
  colour = 'blue';
end;

[U, D] = eig( C );
A = U * sqrtm( D );
b = m;

nTheta = 8;
nPhi   = 8;

nSamples = nTheta * nPhi;

sample = 1;
for t = 1 : nTheta;
  for p = 1 : nPhi;
    theta = 2 * pi * (t / nTheta);
    phi   = pi * ((p / nPhi) - 0.5);
    
    x(1, sample) = cos( phi ) * cos( theta );
    x(2, sample) = cos( phi ) * sin( theta );
    x(3, sample) = sin( phi );
    sample = sample + 1;
  end;
end;

y = A * x + b * ones( 1, nSamples );

for sample = 1:nSamples
  y0(:, 1) = y(:, sample);
  y1(:, 1) = y(:, mod(sample, nSamples) + 1);
  y2(:, 1) = y(:, mod(sample + nPhi - 1, nSamples) + 1);
  
  l1 = [y0'; y1'];
  l2 = [y0'; y2'];
  
if (mod(sample, nPhi) ~= 0)
  line(l1(:, 1), l1(:, 2), l1(:, 3), 'color', colour );
end;
line(l2(:, 1), l2(:, 2), l2(:, 3), 'color', colour );

end;
  








