function y=custom_matvec(matvecdata, transpose)
% y=custom_matvec(matvecdata)
% Custom matrix-vector product function.
%
% INPUT: 
% matvecdata - struct
%            Contains type, specific matrix(~ces), vector(s)
% transpose - if given as "t", uses the transposed matrix
%
% OUTPUT:
% y - vector resulting from the matrix-vector multiplication 

if matvecdata.type == "matrix"
    if nargin == 2 && transpose == "t"
        y=matvecdata.A'*matvecdata.b;
    else
        y=matvecdata.A*matvecdata.b;
    end
    return     
elseif matvecdata.type == "schur_comp"
    if nargin == 2 && transpose == "t"
        y=matvecdata.A'*(matvecdata.M'\(matvecdata.A*matvecdata.b)); 
    else
        y=matvecdata.A'*(matvecdata.M\(matvecdata.A*matvecdata.b)); 
    end
    return  
elseif matvecdata.type == "MI_right_prec_gmres"
%   assumes a right preconditioner blkdiag(M,I)
    [m,~]=size(matvecdata.A);
    b1=matvecdata.b(1:m);
    b2=matvecdata.b(m+1:end);
    y=[b1 + matvecdata.A*b2; ...
        matvecdata.A'*(matvecdata.M\b1)];    return     
elseif matvecdata.type == "saddle"
    [m,~]=size(matvecdata.A);
    b1=matvecdata.b(1:m);
    b2=matvecdata.b(m+1:end);
    y=[matvecdata.M*b1 + matvecdata.A*b2; ...
        matvecdata.A'*b1];    return     
else 
    error("Unknown type. Choose among; matrix, schur_comp, " + ...
        "MI_right_prec_gmres, saddle.")
end


