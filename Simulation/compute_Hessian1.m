%   Compute the Hessian of a real-valued function numerically
%   This is a translation of the Gauss command, hessp(fun,x0), considering only 
%   the real arguments.
%   f: real-valued function (1 by 1)
%   x0: k by 1, real vector
%   varargin: various passing arguments
%   H: k by k, Hessian of f at x0, symmetric matrix


function hf = compute_Hessian1(f,x0,eps,varargin)


% initializations

l_x0=length(x0); % length of x0;

for i=1:l_x0
    x1 = x0;
    x1(i) = x0(i) - eps ;
    df1 = NumJacob(f, x1,varargin{:});
    
    x2 = x0;
    x2(i) = x0(i) + eps ;
    df2 = NumJacob(f, x2,varargin{:});
    
    d2f = (df2-df1) / (2*eps );
    
    hf(i,:) = d2f';
end
end