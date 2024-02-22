function [h,w]=ortho(w,V,method)
% ORTHO:  Orthogonalizes a vector w against a set of vectors V according to
% method.
%
% [h,w]=ortho(w,V,method)
%
% INPUT:
% w - vector to be orthogonalized;
% V - set of vectors to be orthogonalized against;
% method - orthogonalization method:
%           0 - Gram-Schmidt;
%           >=1 - modified Gram-Schmidt, repeated this many times;
%
% OUTPUT:
% h - vector of orthogonalization coefficients;
% w - orthogonalized vector;

if method == 0
    % Gram-Schmidt ortho
    h=V'*w;
    w=w-V*h;
elseif method>0
    colN=size(V,2);
    % (repeated) modified Gram-Schmidt ortho
    h=zeros(colN,1);
    for k=1:method
        for j=1:colN
            v=V(:,j);
            hj=w'*v;
            w=w-hj*v;
            h(j)=h(j)+hj;
        end
    end
else
    error(['Orthogonalization: 0: Gram-Schmidt;', ...
        '>0: iterative modified Gram-Schmidt, run this many times.'])    
end