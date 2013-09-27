%%
%  We study the influence of rho in the flow defined by the Finsler
%  gradient associated with the penalty for piecewise similarity motions.
%  We  perform such an evolution in the case of a rod with different values of rho. 




addpath('toolbox/');

load_mosek();


n = 128*10;

%%
% Helper functions.

% gaussian smoothing
t = [0:n/2, -n/2+1:-1]';
normalize = @(x)x/sum(x);
gauss = @(sigma)normalize( exp( -t.^2/(2*sigma^2) ) );
smooth = @(x,sigma)ifft( fft(x).*fft(gauss(sigma)) );
% project in [0,1]^2
projcurve = @(gamma)clamp(real(gamma),0,1) + 1i*clamp(imag(gamma),0,1);
% recenter
sigmas = 10;  c = (1+1i)/2; rho = .7;
recenter = @(gamma)(gamma-c)*rho + c;

%%
% Load a initial curve

name = 'rod';

% Initial curve

gamma0 = smooth( load_curve(name,n), sigmas );

%%
% Parameter for the energy

% Kernel
options.sigma= .8;
options.delta= .3;


%%
% Helper. 

lw = 2;

cplot = @(gamma,c)plot(real(gamma([1:end 1])), imag(gamma([1:end 1])), 'color', c, 'LineWidth', lw);

%%
% Finsler gradient parameters.

% type of constraint

options.constraints = 2;

% for linprog
options.verbose = 0;
options.linprog_niter = 100;
options.linprog_tol = 1e-12;


%%
% Descent parameters.

% iterations #
niter = 4;
% Display options
cm = jet(niter); % color map
% parameters for linesearch
tau_max = .2; % initial maximum tau
niter_gsec = 5; % # iterations of linesearch


%%
% Gradient

Force = @(gamma)  - 5*real(gamma) - 1000i*(real(gamma)-1/2).^2;

 % descent step
tau = .04;

% number of figures
count=0; 

% constraints parameters
options.lambda=0;
RHO = [0; .3;.7;.9];

%%
% Comparison with different rho

for j=1:4
count=count+1;
gamma = gamma0;
options.rho = RHO(j);  
rho = RHO(j);
figure(count); clf;
for i=1:niter
[op,vec,mat] = load_rigidification_operators(gamma);
grad = perform_rigidification_linprog(gamma, Force(gamma), rho, options);

grad = grad/max(abs(grad)); 

gamma = gamma - tau*(grad);

 figure(count); 
  if mod(i-1,1) == 0
   plot(gamma0); hold on;
   cplot(gamma, cm(i,:)); hold on ;
   axis equal; axis off;
   drawnow; 
  end

end


fig = figure(count);
namenum = [name num2str(count)];
namefile = ['/Users/nardi/Dropbox/banach-optim/NumericalTests/results/compare_rho/' namenum '-rho.eps'];
saveas(fig, namefile);
print('-depsc', namefile);

end

