%%
% Test for Sobolev and Finsler gradient descent
%
% options.gradtype='sobolev' -> Sobolev evolution
% options.gradtype='finsler' -> Finsler evolution 
%
% In the case of the Finsler evolution the parameter option.constraints 
% allows one to choose between the piecewise rigid Finsler gradient 
% (options.constraints=1) and the piecewise similarity Finsler gradient 
% (options.conatrints=2).
%
% We display the evolution, every step of the evolution separately, and
% the graph of the enrgy. 

addpath('toolbox/');
addpath('data/');

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
sigmas = 1;  c = (1+1i)/2; rho = .7;
recenter = @(gamma)(gamma-c)*rho + c;
% evolution
lw = 2;
cplot = @(gamma,c)plot(real(gamma([1:end 1])), imag(gamma([1:end 1])), 'color', c, 'LineWidth', lw);

%%
% Options display

options.lw = 2; % linewidth
options.ncol = 60; % number of colors used
options.nrepeat = 10; % number of times the color patterns repeats
options.mode = 'regular';
options.mode = 'random';

%%
% Load a initial curve

name  = 'harm';

names = {'harm1' 'harm2'};
             
gamma0 = smooth( load_curve(names{1}, n), sigmas );
gamma1 = smooth( load_curve(names{2}, n), sigmas );

%%
% Parameters 

% Kernel

options.sigma=.8;
options.delta=.04;

% Constraints

options.rho=.85;
options.lambda = 2000;

% Energy

Force = @(gamma)compute_rkhs_energy(gamma,gamma1,options);

%%
% Type of descent

%options.gradtype = 'sobolev';
options.gradtype = 'finsler';

niter = 36; % iterations

%%
% Draw the curve.

figure(1); clf; hold on;

plot(gamma1, 'k--' , 'LineWidth', 1.5); hold on;
plot_curve_colored(gamma0, options); hold on;
axis equal; axis off;
drawnow;
fig1 = figure(1);



%%
% Helper. 

lw = 2;

cplot = @(gamma,c)plot(real(gamma([1:end 1])), imag(gamma([1:end 1])), 'color', c, 'LineWidth', lw);

%%
% Finsler gradient parameters.

% type of constraint

options.constraints = 2;% 1= piecewise rigid motion / 2 = piecewise similarity motion

% for linprog
options.verbose = 0;
options.linprog_niter = 100;
options.linprog_tol = 1e-12;

%%
% Sobolev gradient parameter
sigmaS= 25;
options.sigmaS=sigmaS; 

%%
% Descent parameters.

% Display options
kdisp = 1; % frequency of display.
cm = jet(niter); % color map
% parameters for linesearch
tau_max = .2; % initial maximum tau
niter_gsec = 5; % # iterations of linesearch

%%
% Gradient descent

namegrad = [name options.gradtype] ;

namefile = [ '/Users/nardi/Dropbox/banach-optim/NumericalTests/results/harm/' namegrad '-initial.eps'];
saveas(fig1,namefile);
print('-depsc', namefile);

energy=zeros(niter,1);
fig2 = figure(2);
fig3 = figure(3);
fig4 = figure(4);
gamma = gamma0;
Energy = [];
figure(3); clf; hold on;
figure(4); clf;
count =1;

for i=1:niter   
    
    i
   
    % compute the L2 gradient
    [Energy(i), u0] = Force(gamma);
    % modify it to obtain the descent direction
    switch options.gradtype
        case 'sobolev'
            u = smooth(u0, sigmaS);
        case 'finsler'            
            u = perform_rigidification_linprog( gamma, u0, rho, options );
            options.u_init = u; % for next iteration  
            
           [op,vec,mat] = load_rigidification_operators(gamma);
         
    end
    
    u = u/max(abs(u));
    
    % linesearch for the step size
    
    func = @(tau)compute_rkhs_energy(gamma - tau*u, gamma1, options);
    tau = golden_section( func, 0, tau_max, niter_gsec );
    
    if 0
        % display for debug
        tlist = linspace(0,tau_max,10) ;
        E = arrayfun( func, tlist );
        clf; hold on;
        plot(tlist, E); 
        plot(tau, func(tau), 'r*' );
    end    
    if tau==tau_max
        warning('You can increase the step size max');
    end
    
    tau_max = tau*2;
    
    % Descent
    gamma = gamma - tau * u;
    
 
   % plot energy 
    [E,~] = compute_rkhs_energy(gamma,gamma1,options);
    energy(i) = E;
    figure(2); clf;
    axis equal; axis on;
    plot(energy);
    drawnow;
 
   
    % plot evolution 
    
     figure(3); hold on;
     cplot(gamma, cm(i,:)); hold on ;
     axis equal; axis off;
     drawnow; 
     
     figure(4); clf;
    
    % plot iterations   
       
      hiter = figure(count+4); 
      figure(count +4); clf;
      plot(gamma1,'k--', 'LineWidth', 1.5); hold on;
      plot_curve_colored(gamma, options); hold on;
      axis equal; axis off;
      drawnow;
      
      iter = num2str(count);
      nameiter = [namegrad  iter] ;

      namefile = ['/Users/nardi/Dropbox/banach-optim/NumericalTests/results/harm/' nameiter '.eps'];
      saveas(hiter, namefile); 
      print('-depsc', namefile);
      
      count=count+1;

end

fig3 = figure(3);
namefile = ['/Users/nardi/Dropbox/banach-optim/NumericalTests/results/harm/' namegrad '-evolution.eps'];
saveas(fig3, namefile);
print('-depsc', namefile);

fig2 = figure(2);
namefile = ['/Users/nardi/Dropbox/banach-optim/NumericalTests/results/harm/' namegrad '-energy.eps'];
saveas(fig2, namefile);
print('-depsc', namefile);


figure(4); clf;

plot(gamma1, 'k--', 'LineWidth', 1.5);
hold on;
plot_curve_colored(gamma,options);
hold on;
axis equal; axis off;
drawnow;
namefile = ['/Users/nardi/Dropbox/banach-optim/NumericalTests/results/harm/' namegrad '-final.eps'];
saveas(fig4, namefile);
print('-depsc', namefile);








