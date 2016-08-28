function [u,err] = perform_rigidification_linprog(gamma, u0, rho, options)

%  perform_rigidification_linprog - Finsler gradient using linear programming
%
%   u = perform_rigidification_linprog(gamma, u0, rho, options);
%   
%   1) Piecewise rigid motions (options.constraints=1)
%
%   Solves the problem
%
%       min_u |H(u)|_1 s.t. |<u-u0,N>|<=rho*|<u0,N>| and L(u)=0.
%
%   A large  rho in [0,1] corresponds to a strong rigidification. 
%   
%  2) Piecewise similarity motions (options.constraints=2)
%
%   Solves the problem
%
%       min_u |B(u)|_1 s.t. |<u-u0,N>|<=rho*|<u0,N>| and |L|_1<=lambda
%
%   A large  rho in [0,1] corresponds to a strong rigidification.
%   A large positif lambda authorizes a large scaling
%
%   Copyright (c) 2012 Gabriel Peyre



%%
% Models: rigid (constraints=1) or similarity (constraints=2) 

constraints = getoptions(options, 'constraints', 'inf');

%%
% The data to be rigidified

[op,vec,mat] = load_rigidification_operators(gamma);

x0 = op.Vs(u0);


%%
% Mosek form
%   min  <c,weight>    s.t.
%   Amin<=A*w<=Amax,
%   Wmin<=w<=Wmax


options.null = 0;
n = length(gamma);

switch constraints
        
    case {'1' 1 'l1'} % Piecewise rigid motions
        
        % value of the constraint
        epsilon = norm(real(x0(:)), 1)*rho;
 
        
        % Solve the linear program w=[gamma;t;s]
        %       where t has size nH and s has size n
        % min <t,weight>  s.t.
        %   0 <= L*gamma   <=0          (C1)        n
        %        H*gamma-t <=0,         (C2)        n
        %   0 <= H*gamma+t              (C3)        n
        %        <gamma,N>-s <= x0      (C4)        n
        %   x0<= <gamma,N>+s            (C5)        n
        %           <s,1>    <= epsilon (C6)        1
        
        A = [   mat.L,  sparse(n,n), sparse(n,n); ...         % (C1)
            mat.H, - speye(n,n), sparse(n,n); ...             % (C2)
            mat.H, + speye(n,n), sparse(n,n); ...             % (C3)
            mat.dp(vec.N), sparse(n,n), -speye(n,n); ...      % (C4)
            mat.dp(vec.N), sparse(n,n), +speye(n,n); ...      % (C5)
            zeros(1,2*n),    zeros(1,n),   ones(1,n);   ...   % (C6)
            ];
        Amin = [    zeros(n,1); ...     % (C1)
            -Inf + zeros(n,1); ...      % (C2)
            zeros(n,1);   ...           % (C3)
            -Inf + zeros(n,1);   ...    % (C4)
            x0; ...                     % (C5)
            -Inf  ...                   % (C6)
            ];
        Amax = [    zeros(n,1); ...     % (C1)
            zeros(n,1); ...             % (C2)
            +Inf + zeros(n,1); ...      % (C3)
            x0; ...                     % (C4)
            +Inf + zeros(n,1); ...      % (C5)
            epsilon  ...                % (C6)
            ];
        Xmin = [];
        Xmax = [];
        C = [zeros(2*n,1); ones(n,1); zeros(n,1)];
        
        
        case {'1' 2 'l1'} % Piecewise similarity motions
            
        % value of the constraint
        epsilon = norm(real(x0(:)), 1)*rho;
        
        % Solve the linear program w=[gamma;t;z;s]
        %       where t has size nH and s has size n z has size nL
        % min <t,1>    s.t.
        %        
        %        B*gamma -t  <=0        (C1)        2*n
        %   0 <= B*gamma+t              (C2)        2*n
        %        L*gamma-z <= 0,        (C3)        n
        %   0 <= L*gamma+z              (C4)        n
        %    0<= <z,1> <= lambda        (C5)        1
        %        <gamma,N>-s <= x0      (C6)        n
        %   x0<= <gamma,N>+s            (C7)        n
        %           <s,1>    <= epsilon (C8)        1
        
        A = [
            mat.B, - speye(2*n,2*n), sparse(2*n,n),  sparse(2*n,n) ;... % (C1)
            mat.B, + speye(2*n,2*n), sparse(2*n,n),  sparse(2*n,n) ;... % (C2)
            mat.L,   sparse(n,2*n),   - speye(n,n), sparse(n,n); ...    % (C3)
            mat.L,   sparse(n,2*n),   + speye(n,n), sparse(n,n); ...    % (C4)
            zeros(1,2*n),  zeros(1,2*n),   ones(1,n), zeros(1,n); ...   %( C5)
            mat.dp(vec.N), sparse(n,2*n), sparse(n,n), -speye(n,n) ; ...% (C6)
            mat.dp(vec.N), sparse(n,2*n), sparse(n,n), +speye(n,n); ... % (C7)
            zeros(1,2*n),  zeros(1,2*n),   zeros(1,n), ones(1,n);   ... % (C8)
            ];
        Amin = [ -Inf +  zeros(2*n,1); ... % (C1)
             zeros(2*n,1); ...             % (C2)
            -Inf + zeros(n,1); ...         % (C3)
            zeros(n,1);   ...              % (C4)
            0;...                          % (C5)
            -Inf + zeros(n,1);   ...       % (C6)
            x0; ...                        % (C7)
            -Inf  ...                      % (C8)
            ];
        Amax = [  zeros(2*n,1); ...      % (C1)
            +Inf + zeros(2*n,1); ...     % (C2)
            zeros(n,1); ...              % (C3)
            +Inf + zeros(n,1); ...       % (C4)
            options.lambda; ...          % (C5)
            x0; ...                      % (C6)
            +Inf + zeros(n,1); ...       % (C7)
            epsilon  ...                 % (C8)
            ];
        Xmin = [];
        Xmax = [];
        C = [zeros(2*n,1); ones(2*n,1); zeros(n,1); zeros(n,1)];
        
        
    otherwise
        error('Unknown constraints type');
end


%%
% Setup Mosek variables

prob.c = C;
prob.a = A;
prob.blc = Amin;
prob.buc = Amax;
prob.blx = Xmin;
prob.bux = Xmax;

%%
% Set parameters

param = [];
% max numer of iterations
param.MSK_IPAR_INTPNT_MAX_ITERATIONS = getoptions(options, 'linprog_niter', 100);
% tolerance, primal
param.MSK_DPAR_INTPNT_TOL_PFEAS = getoptions(options, 'linprog_tol', 1e-12);
param.MSK_DPAR_INTPNT_TOL_REL_GAP = getoptions(options, 'linprog_tol', 1e-12);
% verbosity level, 0=nothing is echoed, 3=all is echoed
verb = getoptions(options, 'verbose', 0);

% Perform the optimization.
[r,res] = mosekopt(['echo(' num2str(verb) ') minimize info'], prob, param);
if r~=0
    warning(['Mosek problem: ' res.rcodestr]);
end
err.niter = res.info.MSK_IINF_INTPNT_ITER;
sol = res.sol.itr;
w   = sol.xx;
u = w(1:n) + 1i*w(n+1:2*n);


