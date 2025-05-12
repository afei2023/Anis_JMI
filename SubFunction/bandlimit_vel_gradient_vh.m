function grad_bl=bandlimit_vel_gradient_vh(grad,fmax,vel,dx,dz)
% fancy smoother ... gaussian moving average 
% grad is total energy at each model gridpoint

nx=size(grad,1);
nz=size(grad,2);
dkz=2*pi/(nz*dz); % radial domain
dkx=2*pi/(nx*dx); % radioal domain ... like omega
nkz=ceil(nz/2); % round up nz/ 2 ... to calc half .. and conj ..4 optimize
nkx=ceil(nx/2); % round up nx/2  

kmax=2*pi*fmax/min(min(vel)); % maxiumum ko

% define vertical filter
Filterz=zeros(1,nz);
for ikz=1:nkz
    kz=(ikz-1)*dkz;
    Filterz(ikz)=exp(-kz^2/(0.2.*kmax^2));
end
Filterz(nz-nkz+2:nz)=conj(Filterz(nkz:-1:2)); % other half
filterz=real(fft(Filterz))/nz;
% filterz(nkz:nz)=0;
filterz=circshift(filterz,[1 -nkz]); % to make filter at middel

% define horizontal filter ... same thing for x
Filterx=zeros(1,nx);
for ikx=1:nkx
    kx=(ikx-1)*dkx;
    Filterx(ikx)=exp(-kx^2/(0.2.*kmax^2));
end
Filterx(nx-nkx+2:nx)=conj(Filterx(nkx:-1:2));
filterx=real(fft(Filterx))/nx;
filterx=circshift(filterx,[1 -nkx]);

% extend model
gradtmp=zeros(nx+2*nkx,nz+2*nkz);
gradtmp(nkx+1:nkx+nx,nkz+1:nkz+nz)=grad; % put energy matix in the middel 

% gradtmp(1:nkx,nkz+1:nkz+nz)=repmat(grad(1,:),[nkx 1]);
% gradtmp(end-nkx+1:end,nkz+1:nkz+nz)=repmat(grad(nx,:),[nkx 1]);
% gradtmp(nkx+1:nkx+nx,1:nkz)=repmat(grad(:,1),[1 nkz]);
% gradtmp(nkx+1:nkx+nx,end-nkz+1:end)=repmat(grad(:,nz),[1 nkz]);

% gradtmp(1:nkx,1:nkz)=grad(1,1);
% gradtmp(end-nkx+1:end,1:nkz)=grad(nx,1);
% gradtmp(1:nkx,end-nkz+1:end)=grad(1,nz);
% gradtmp(end-nkx+1:end,end-nkz+1:end)=grad(nx,nz);


grad_bltmp1=zeros(nx+2*nkx,2*nz+2*nkz-1);
% smooth along x
for ix=1:nx+2*nkx
    grad_bltmp1(ix,:)=conv(gradtmp(ix,:),filterz(nz:-1:1)); % conv for smoother
end

grad_bltmp2=zeros(2*nx+2*nkx-1,2*nz+2*nkz-1);
% smooth along z
for iz=1:2*nz+2*nkz-1
    grad_bltmp2(:,iz)=conv(grad_bltmp1(:,iz),filterx(nx:-1:1));
end
grad_bl=grad_bltmp2(2*nkx:2*nkx+nx-1,2*nkz:2*nkz+nz-1);
