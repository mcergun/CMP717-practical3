% PlotColourSpace.m

function PlotColourSpace_vs2(FData, BData, pF, mF, CF, pB, mB, CB);

nSamples = 200;

figure(1001)
close(1001)
figure(1000)
close(1000)

% store the axis 
save_axis = gca;

[hi1 hi2] = size(FData);
[hi3 hi4] = size(BData);
if ((hi1>0) & (hi3>0))

% show figures 
figure(1000);
% title('Data in Colour Space ----  x(red) foreground; *(blue) background');
grid on;
hold on;
axis equal;
axis vis3d;

[nDimensions, nFData] = size(FData);
[nDimensions, nBData] = size(BData);

% get number of Samples

if (nFData < nSamples)
  nSamples = nFData-1;
end
if (nBData < nSamples)
  nSamples = nBData-1;
end

disp(['Number of Samples: ' num2str(nSamples)])
%nSamples = 1;

FPerm = randperm(nFData);
FSamples = FPerm(1:nSamples);
for i = 1:nSamples
  datai = FData(:, FSamples(i));
  plot3(datai(1), datai(2), datai(3), 'x', 'color', datai);
end;

BPerm = randperm(nBData);
BSamples = BPerm(1:nSamples);
for i = 1:nSamples
  datai = BData(:, BSamples(i));
  plot3(datai(1), datai(2), datai(3), '.', 'color', datai);
end;  

% load a good 3dviewer
view3d rot;

end

figure(1001);
title('Gaussian Mixture Model - blue(x) foreground; red(*) background');
grid on;
box on;
hold on;
axis equal;
axis vis3d;

% Info
nF = size(pF,2); % GMMOptions.nComponents;
nB = size(pB,2); % GMMOptions.nComponents;

for i = 1:nF
    DrawEllipse3D(mF{i}, CF{i}, [0 0 1]'); % mF{i}); 
end

for i = 1:nB
    DrawEllipse3D(mB{i}, CB{i}, [1 0 0]'); % mB{i}); 
end

% for i = 1:nB
%  DrawEllipse3D(mM{nalpha}{i}, CM{nalpha}{i}, mM{nalpha}{i});
% end;

% load a good 3dviewer
view3d rot;

% write it
% print -depsc2 colour;
% disp('I wrote the images in the current directory as: colour.eps');    
% pause

% give back the axis
axes(save_axis); 
