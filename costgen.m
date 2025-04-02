function [H,h,const]=costgen(T,S,Q,R,P,dim,x0,xr,ur)

Qbar = blkdiag(kron(eye(dim.N - 1),Q),P); 
Rbar = kron(eye(dim.N),R);
H = S'*Qbar*S+Rbar;

hx0 = S'*Qbar*T*x0;
hxr = -S'*Qbar*kron(ones(dim.N,1),eye(dim.nx))*xr;
hur = -Rbar*kron(ones(dim.N,1),eye(dim.nu))*ur;
h = hx0 + hxr + hur;
const=x0'*T'*Qbar*T*x0;

end