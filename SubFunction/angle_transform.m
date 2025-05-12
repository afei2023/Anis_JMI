%% generate transform matrix for angle-transform

function A=angle_transform(c,alpha,nx,dx,ixmid,om)

x=(0:dx:(nx-1)*dx)-(ixmid-1)*dx;

nom=length(om);
eta=sin(pi*alpha/180);
neta=length(eta);

om3D=repmat(om,[nx 1 neta]); 
om3D=permute(om3D,[3 1 2]);
eta3D=repmat(eta.',[1 nx nom]);
x3D=repmat(x,[neta 1 nom]);

A=exp(-1i*eta3D.*om3D.*x3D./c);
