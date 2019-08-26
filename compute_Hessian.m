
%   Compute the Hessian of a real-valued function numerically
%   This is a translation of the Gauss command, hessp(fun,x0), considering only 
%   the real arguments.
%   f: real-valued function (1 by 1)
%   x0: k by 1, real vector
%   varargin: various passing arguments
%   H: k by k, Hessian of f at x0, symmetric matrix




function H = compute_Hessian(f,x0,eps,varargin)


% initializations

k=length(x0);
hessian=zeros(k,k);
grdd=zeros(k,1);


% Computaion of stepsize (dh)

ax0=abs(x0);
for i=1:k
    if x0(i,1)~=0
       dax0(i,1)=x0(i,1)/ax0(i,1);
    else
       dax0(i,1)=1;
    end
end

dh=eps*max([ax0, (1e-2)*ones(length(x0),1)]')'.*dax0;
xdh=x0+dh;
dh=xdh-x0;  % This increases precision slightly
ee=[]; I=eye(k);
for i=1:k
    ee(:,i)=I(:,i).*dh;
end


% Computation of f0=f(x0)

f0=feval(f,x0,varargin{:});

% Compute forward step

for i=1:k
    grdd(i,1)=feval(f,x0+ee(:,i),varargin{:});
end

% Compute 'double' forward step

for i=1:k
    for j=i:k
        hessian(i,j)=feval(f,x0+(ee(:,i)+ee(:,j)),varargin{:});
        if i~=j
            hessian(j,i)=hessian(i,j);
        end
    end
end

l=ones(1,k); 
grdd=kron(l,grdd); 

H=(((hessian - grdd) - grdd') + f0) ./ (kron(dh,dh'));
end