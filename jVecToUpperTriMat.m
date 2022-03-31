function m=jVecToUpperTriMat(v,size_m)
% converts vector to the upper-triangular part of a matrix 
%
% IN:
%   v: vector
%   size_m: size of matrix (size_m x size_m)
% OUT:
%   m: matrix
% 
% v1.0 Feb 2015 Dimitri Van De Ville
% - initial version

% get indices of upper triangular part (Peter Acklam's trick)
idx = find(triu(ones(size_m), 1));

m=zeros(size_m,size_m);
m(idx)=v(:);
%m=m+m.';
