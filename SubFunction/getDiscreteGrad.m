function [G,Gx] = getDiscreteGrad(n,h)
%n1 is the number of grid points in the vertical (depth) dimension
%n2 is the number of grid points in the lateral dimension
%h1 is the distance between grid points in the vertical (depth) dimension
%(assumed constant)
%h2 is the distance between grid points in the lateral dimension
%(assumed constant)
% dwp are the depth weights. Set to all ones for no depth weighting

dpw=ones(n(1),1);
% dpw_ave = (dpw(1:end-1)+dpw(2:end))/2;
n1=n(1);
n2=n(2);
h1=h(1);
h2=h(2);
%define difference matrix D acting on vectorized model using kron
%I1 = speye(n1); %z
I2 = speye(n2); %x
Gz = spdiags(ones(n1,1)*[-1 1],0:1,n1,n1)/h1;
Gz(end,:)=-Gz(end-1,:);
% Gx = spdiags(ones(n2,1)*[1 -1],-1:0,n2,n2-1)/h2;
Gx = spdiags(ones(n2,1)*[1 -1],-1:0,n2,n2)/h2;

if nargout==1
    G = [kron(I2,spdiags(dpw,0,n1-1,n1-1)*Gz) ;kron(Gx',spdiags(dpw,0,n1,n1))]; 
    %DTD is also minus discrete Laplacian % x and y derivatives in modified Laplacian
else
	G =kron(I2,spdiags(dpw,0,n1,n1)*Gz);		% z derivatives 
    % [G] = insertrow(G,n2,1);
    Gx1 =kron(Gx,spdiags(dpw,0,n1,n1));
	Gx =kron(Gx',spdiags(dpw,0,n1,n1));					% x derivatives 
    Gx(end-n1+1:end,:)=Gx1(end-n1+1:end,:);
    % [Gx] = insertrow(Gx,n1,0);

end


% E = (D<0); %keep track of forward differences for defining vector l1 norm for isotropic TV

end