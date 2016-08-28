%%
%  We compare the evolution by using different gradients: 
%  1) L^2 gradient;
%  2) Finsler gradient using the penalty for piecewise rigid motions: this
%     corresponds to calculate the Finsler gradient associated with the
%     penalty for piecewise similarity motions setting lambda = 0;
%  3) Finsler gradient using the penalty for piecewise similarity motions.
%  


addpath('toolbox/');

load_mosek();

%%
% Number of points

n = 128*10;

%%
% Helper functions.

% gaussian smoothing, only for periodic curves 
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

% Gradient

Force = @(gamma) - 5*real(gamma) - 1000i*(real(gamma)-1/2).^2;

%%
% Helper. 

lw = 2;

cplot = @(gamma,c)plot(real(gamma([1:end 1])), imag(gamma([1:end 1])), 'color', c, 'LineWidth', lw);

%%
% Finsler gradient parameters.

% type of constraint

options.constraints = 2;% 1= piecewise rigid motions / 2 = piecewise similarity motions

% for linprog
options.verbose = 0;
options.linprog_niter = 100;
options.linprog_tol = 1e-12;
 

%%
% Descent parameters.

% iterations 
niter = 16;

% Display options
cm = jet(niter); % color map
% parameters for linesearch
tau_max = .2; % initial maximum tau
niter_gsec = 5; % # iterations of linesearch

%%
%  Gradients comparison

gamma1 = gamma0;
gamma2 = gamma0;
gamma3 = gamma0;

%%
% Descente step

tau = .03;


%%
% Numerical tests

% Gradient: L2
figure(1); clf;
for i=1:niter
   grad = Force(gamma1);
   grad = grad/max(abs(grad)); 
   gamma1 = gamma1 - tau*grad;
   figure(1);
   if mod(i-1,5) == 0
   plot(gamma0); hold on;
   cplot(gamma1, cm(i,:)); hold on ;
   axis equal; axis off;
   drawnow;
   end
end
fig1 = figure(1);
namefile = ['/Users/nardi/Dropbox/banach-optim/NumericalTests/results/compare_gradients/' name '-L2.eps'];
saveas(fig1, namefile);
print('-depsc', namefile);


% Gradient: Finsler lambda =0 rho=.5
figure(2); clf;
options.lambda = 0;
rho = .5;
options.rho = rho ;
for i=1:niter
  grad = perform_rigidification_linprog(gamma2, Force(gamma2), rho, options);
  grad = grad/max(abs(grad)); 
  gamma2 = gamma2 - tau*grad;
  figure(2);
    if mod(i-1,5) == 0
   plot(gamma0); hold on;
   cplot(gamma2, cm(i,:)); hold on ;
   axis equal; axis off;
   drawnow;
    end
end
fig2 = figure(2);
namefile = ['/Users/nardi/Dropbox/banach-optim/NumericalTests/results/compare_gradients/' name '-Finsler1.eps'];
saveas(fig2, namefile);
print('-depsc', namefile);

% Gradient: Finsler weightL >> 1 rho=.5
figure(3); clf;
options.lambda = 2000;
rho = .5;
options.rho = rho ;
for i=1:niter
 grad = perform_rigidification_linprog(gamma3, Force(gamma3), rho, options);
 grad = grad/max(abs(grad)); 
 gamma3 = gamma3 - tau*grad;
 figure(3);
   if mod(i-1,5) == 0
   plot(gamma0); hold on;
   cplot(gamma3, cm(i,:)); hold on ;
   axis equal; axis off;
   drawnow; 
   end
end
fig3 = figure(3);
namefile = ['/Users/nardi/Dropbox/banach-optim/NumericalTests/results/compare_gradients/' name '-Finsler2.eps'];
saveas(fig3, namefile);
print('-depsc', namefile);

