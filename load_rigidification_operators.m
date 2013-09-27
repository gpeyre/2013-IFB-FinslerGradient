
function [op,vec,mat] = load_rigidification_operators(gamma)

% load_rigidification_operators  loads the linear operator associated to
% the calculation of the finsler gradient
%
%   [op,vec,mat] = load_rigidification_operators(gamma);
%
%   vec.N/vec.T is the normal computed by 2nd order differences
%       (one per vertex)
%   vec.T1/vec.N1 is the tangent/normal computed by 1st order differences 
%       (one per edge).
%
%       op.V(x) = real(x)*N + imag(x)*T
%       op.H(u) = D^*( ( D(u) ) . N1 )
%       op.L(u) = D(u) . T1
%
%   where "." is the point-wise inner product (op.dotp)
%   D is the forward derivative (op.DerF)
%   D^* is the backward derivative (op.DerB)
%   
%   op.Vs, op.Hs, op.Ls are the adjoints. 
%
%   mat.H and mat.L contains the (sparse) matrix version, 
%   taking 2*n size real/complex input [real(u);imag(u)]
%
%   Copyright (c) 2012 Gabriel Peyre


%%
% Helpers functions

normalize = @(c)c./max(abs(c),1e-10);
rotate = @(c)1i*c;
% pointwise dotproduct
dotp = @(x,y)real(x.*conj(y));


%% 
% Derivatives operators, tangent and normal vectors.

    DerF    = @(x)x([2:end 1]) - x;        % first order derivative
    DerB    = @(x)x - x([end 1:end-1]);    % -adjoint    
    DerC    = @(x)(x([2:end 1]) - x([end 1:end-1]))/2;
    DerC0   = DerC;
    Der2    = @(x)DerF(DerB(x));   % 2nd order derivative
    Der2S   = @(x)Der2(x);         % adjoint


% tangent and normal at each vertex
T = normalize(DerC0(gamma));
N = rotate(T);
% first order tangent (on edge)
T1 = normalize(DerF(gamma));
N1 = rotate(T1);

%%
% V L H operators
V = @(x)real(x).*N + imag(x).*T;
Vs = @(v)dotp(v,N) + 1i*dotp(v,T);

L  = @(u)dotp(DerF(u),T1);
Ls = @(w)-DerB(w.*T1);
    
H = @(u)DerB(dotp(DerF(u),N1));    
Hs = @(w)DerB( DerF(w).*N1 );

%%
% Export symbols.

op.normalize = normalize;
op.rotate = rotate;
op.dotp = dotp; 
%
op.V    = V;
op.Vs   = Vs;
op.L    = L;
op.Ls   = Ls;
op.H    = H;
op.Hs   = Hs;
%
op.DerF    = DerF;
op.DerB    = DerB;
op.DerC    = DerC;
op.Der2    = Der2;
op.Der2S   = Der2S;
op.DerC0   = DerC0;
%
vec.T   = T;
vec.N   = N;
vec.T1  = T1;
vec.N1  = N1;

% inner product with complex vector w
mat.dp = @(w)[spdiags(real(w),0,length(w),length(w)) spdiags(imag(w),0,length(w),length(w))];

n=size(gamma, 1);
% forward derivative
mat.DerF = -speye(n) + spdiags(ones(n,1),1,n,n);
mat.DerF(n,1)=1;
% backward derivative
mat.DerB = -mat.DerF';
% forward derivative applied to both component
mat.DerFVec = [mat.DerF, sparse(n,n); sparse(n,n), mat.DerF];
% L, H operators
mat.L = mat.dp(vec.T1) * mat.DerFVec;
mat.H = mat.DerB * mat.dp(vec.N1) * mat.DerFVec;
% B operator
mat.HT = mat.DerB * mat.dp(vec.T1) * mat.DerFVec;
mat.B = [mat.HT; mat.H];
        



end

