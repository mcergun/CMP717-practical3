% [im_to] = image_index_const(index,im_to,im_from); 
%
% Method: im_to(index) = im_from; 
%         with im_from = 3x1 vector
% 

function [im_to] = image_index(index,im_to,im_from)

% decompose all
im_to1 = im_to(:,:,1);
im_to2 = im_to(:,:,2);
im_to3 = im_to(:,:,3);

% assigne it 
im_to1(index) = im_from(1);
im_to2(index) = im_from(2);
im_to3(index) = im_from(3);

% compose it 
im_to(:,:,1) = im_to1;
im_to(:,:,2) = im_to2;
im_to(:,:,3) = im_to3;
