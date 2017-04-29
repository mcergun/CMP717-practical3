%------------------------------------------------------------------
% This is a simple example of finding a white object on a dark
% background in an image, where you impose that spatial coherence
%
% Currently it find the global optimujm of the enrgy:
%
% E = Sum_(all pixel i) Unary term(i,li) + Sum_(all pairs of pixels i,j) 
%     [li diff lambda1 lj] ( lambda1 + lambda2 * exp - beta ||x_i-x_j||^2 )
%
% where li is the label of pixel i, and x_i is the colour (in 0,1 range) at pixel i, and [a] is 1 if arguemnt a is true 
%
%------------------------------------------------------------------

function GCPreparation();

% Settings
Name = 'flower.jpg'; % Square.png
Parameters.lambda1 = 1.0; % weight of the Ising Prior
Parameters.neighborSystem = 1; % 0-4Neig.; 1-8Neig. - where the diagonal terms are weighted with a factor of 1/sqrt(2)
Parameters.lambda2 = 10.0; % weight of the Edge Dependent Term
output = 0;
numLabels = 2; % number of labels

% load the image
image = double(imread(Name))/255;
[M,N,S] = size(image);

% prepare the Unary Potentials
UnaryPotentials = zeros(M,N,numLabels);

% compute the Unary term Fgd is white, Bkg is dark
hm1 = zeros(M,N,3);
hm2 = abs(image-hm1);
UnaryPotentials(:,:,1) = (sum(hm2,3))/3;
hm1 = ones(M,N,3);
hm2 = abs(image-hm1);
UnaryPotentials(:,:,2) = (sum(hm2,3))/3;

% compute beta
count = 0;
beta = 0.0;
hm1 = image;
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
tic
disp('Call GraphCut ... ');
LabelOut = GraphCut(UnaryPotentials,Parameters.lambda1,Parameters.lambda2,0,Parameters.neighborSystem,image,beta);
toc


% show the results
figure(1);
imshow(image);
figure(2);
imshow(LabelOut);
figure(3);
hm2 = zeros(M,N,3);
hm2(:,:,1) = LabelOut;
hm2(:,:,2) = LabelOut;
hm2(:,:,3) = LabelOut;
hm1 = 0.5*image + 0.5*hm2;
imshow(hm1);

% % TEST
% count = 1;
% for hi1 = 0:1/255:1
%     
%     d1 = (0-hi1)^2;
%     d1 = d1+ (0-hi1)^2;
%     d1 = d1+ (0-hi1)^2;
%     
%     d2 = 1+50*exp(-beta*d1);
%     k1(count) = d2;
%     
%     d3 = 1.0/(0.02+2.0*sqrt(d1));
%     k2(count) = d3;
%     
%     count = count+1;
% end
% 
% figure(10)
% clf;
% plot(k1)
% hold on;
% plot(k2)
