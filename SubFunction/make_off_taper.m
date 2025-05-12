function Ptaper=make_off_taper(noff,sxloc,nx,nxtap)

nsrc=length(sxloc);
nx2=nx+2*nxtap;

Ptaper=zeros(nx2,nsrc);
for ix=1:nsrc
    ixsrc=sxloc(ix)+nxtap;
    x1=ixsrc-noff;
    if x1 > 0
        Ptaper(x1:ixsrc,ix)=1;
    else
        Ptaper(1:ixsrc,ix)=1;
    end
    
    x2=ixsrc+noff;
    if x2 < nx2
        Ptaper(ixsrc:x2,ix)=1;
    else
        Ptaper(ixsrc:nx2,ix)=1;
    end
end

