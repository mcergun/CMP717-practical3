% example to compute likelihoods

nameIm = '../data/starfish.bmp';
nameUser = '../data/starfishUser.bmp';

%-----------------------------------
% Settings

addpath 'GMM';
addpath 'Optimization';
addpath 'useful';

% general
visualizeIntermediateResults = 1;
visualizeFinalResults = 1;

%--------------------------
% START: Modify the lines below

% some MRF settings
Parameters.neighborSystem = 1; % CHOOSE 0 or 1; (0-4Neighborhood; 1-8Neighborhood)
Parameters.lambda1 = 0.0; % CHOOSE > 0. weight of the Ising Prior
Parameters.lambda2 = 0.0; % weight of the Edge Dependent Term

% END: Modify the lines above
%--------------------------

% some GMM settings
optionsGMM.nComponents = 10; 
optionsGMM.iterations = 2;
optionsGMM.init = 0; 
optionsGMM.method = 'kmeans';                      				    
optionsGMM.covarType = 'full';
optionsGMM.min_var = (1/255)^2.0;

extra.ValueInf = 100000;

%------------------------------------


% load data
im = double(imread(nameIm))/255;
[M,N,S] = size(im);
imr = im(:,:,1);
img = im(:,:,2);
imb = im(:,:,3);
input = imread(nameUser);
inputr = input(:,:,1);
inputg = input(:,:,2);
inputb = input(:,:,3);


%------------------------------
% LIKELIHOOD

% get GMMs
disp('compute likelihood from brushes ... ');
indexF = find(inputr==255 & inputg==0 & inputb==0);
indexB = find(inputr==0 & inputg==0 & inputb==255);

dataF = zeros(3,size(indexF,1));
dataF(1,:) = imr(indexF);
dataF(2,:) = img(indexF);
dataF(3,:) = imb(indexF);

dataB = zeros(3,size(indexB,1));
dataB(1,:) = imr(indexB);
dataB(2,:) = img(indexB);
dataB(3,:) = imb(indexB);

% unused values
p =0; m =0; C = 0; sigma = 0;
init = 1; % initialize always
[pF, mF, CF] = ComputeGMM_cr( dataF, optionsGMM, init, p, m, C, sigma);
[pB, mB, CB] = ComputeGMM_cr( dataB, optionsGMM, init, p, m, C, sigma);


% visualize GMMs
if (visualizeIntermediateResults)
    PlotColourSpace_vs2(dataF, dataB, pF, mF, CF, pB, mB, CB);
end

% now evaluate the GMM
data = [imr(:)';img(:)';imb(:)']; % in vector form
fLikeVec = EvGMM_new(data, pF, mF, CF);
bLikeVec = EvGMM_new(data, pB, mB, CB);

fLike = reshape(fLikeVec,[M,N]);
bLike = reshape(bLikeVec,[M,N]);

% visualize the likelihood
if (visualizeIntermediateResults)
    figure(100)
    clf;
    hm1 = fLike - bLike; 
    imagesc(hm1)
    colorbar;
    colormap(jet)
    title('-log likelihood Fgd-Bgd (ratio)');

    figure(101)
    clf;
    imagesc(hm1<0);
    title('optimal result with likelihood only');

end


%------------------------------
% set up the energy

% set user hard user constraints
bLike(indexF) = extra.ValueInf;
bLike(indexB) = -extra.ValueInf;
fLike(indexB) = extra.ValueInf;
fLike(indexF) = -extra.ValueInf;

% prepare the Unary Potentials
UnaryPotentials = zeros(M,N,2);
UnaryPotentials(:,:,1) = bLike;
UnaryPotentials(:,:,2) = fLike;

% compute beta
count = 0;
beta = 0.0;
hm1 = im;
hm2 = (hm1(2:M,:,:)-hm1(1:(M-1),:,:)).^2;
[c1 c2 c3] = size(hm2);
count = count + (c1*c2*c3);
beta = beta + sum(sum(sum(hm2)));
hm2 = (hm1(:,2:N,:)-hm1(:,1:(N-1),:)).^2;
[c1 c2 c3] = size(hm2);
count = count + (c1*c2*c3);
beta = beta + sum(sum(sum(hm2)));
if (Parameters.neighborSystem)
    hm2 = (hm1(2:M,2:N,:)-hm1(1:(M-1),1:(N-1),:)).^2;
    [c1 c2 c3] = size(hm2);
    count = count + (c1*c2*c3);
    beta = beta + sum(sum(sum(hm2)));
    hm2 = (hm1(1:(M-1),2:N,:)-hm1(2:M,1:(N-1),:)).^2;
    [c1 c2 c3] = size(hm2);
    count = count + (c1*c2*c3);
    beta = beta + sum(sum(sum(hm2)));
end

% final beta
beta_inv = 2.0 * (beta/count);
if (beta_inv == 0)
    beta_inv = epsilon;
end
beta = 1.0 / beta_inv;
disp(['The Edge Beta is: ',num2str(beta)]);

% Call Graph Cut
disp('Call GraphCut ... ');
LabelOut = GraphCut(UnaryPotentials,Parameters.lambda1,Parameters.lambda2,0,Parameters.neighborSystem,im,beta);

% show the final results
if (visualizeFinalResults)
    
    figure;
    clf;
    bkgIm = ones(M,N,3);
    bkgIm(:,:,3) = 1;
    index = find(LabelOut==1);
    bkgIm = image_index(index,bkgIm,im);
    imagesc(bkgIm);

    title('final segmentation result');
    
end

  












