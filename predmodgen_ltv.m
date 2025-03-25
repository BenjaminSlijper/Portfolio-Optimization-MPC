function [P,S]=predmodgen(sys, k0 dim)

%Prediction matrices generation
%This function computes the prediction matrices to be used in the
%optimization problem

%Prediction matrix from initial state
P=zeros(dim.ny*(dim.N),dim.nx);
[A,B,...]=sys.DataFcn(k0)
cumA = A;
P(1:dim.ny,:)=C*A;
for k=1:dim.N-1
    [A,B,...]=sys.DataFcn(k0+k)
    P(k*dim.ny+1:(k+1)*dim.ny,:)=C*cumA*A;
    cumA=cumA*A;
end

%Prediction matrix from input
S=zeros(dim.ny*(dim.N),dim.nu*(dim.N));
for k=1:dim.N-1
    for i=0:k-1
        S(k*dim.ny+1:(k+1)*dim.ny,i*dim.nu+1:(i+1)*dim.nu)=C*LTI.A^(k-1-i)*B;
    end
end

