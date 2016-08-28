function plot_curve_colored(gamma, options)

% plot_curve_colored - plot a periodic curve with colors
%
%   plot_curve_colored(gamma, options);
%
%   If options.nrepeat>0 and options.p>0: uses random piecewise constant colors 
%
%   If options.err not empty: uses err in [0,1] to index the blue/red
%   colors.
%
%   If options.colors is provided: uses these colors.
%
%   Copyright (c) 2013 Gabriel Peyre


n = length(gamma);

lw  = getoptions(options, 'lw', 2);
p = getoptions(options, 'ncol', -1);
t  = getoptions( options, 't', (0:n-1)/n );
nrepeat = getoptions(options, 'nrepeat', -1);
mode = getoptions(options, 'mode', 'regular');
err = getoptions(options, 'err', []);
Col = getoptions(options, 'colors', []);

if nrepeat>0 && p>0
    %%% MODE RANDOM COLORS %%%
    
    % colors
    r = 512;
    u = linspace(0,1,r)';
    C = hsv2rgb( [u u*0+1 u*0+1] );
    j = 1;
    
    jump = r*nrepeat / p;
    
    % set the random number to fixed initialization
    rand('state',0);
    
    switch mode
        case 'random'
            J = 1+floor(rand(p,1)*(r-1));
        case 'semi-random'
            %%%  Random mode %%%
            J = j;
            for i=1:p
                j = j + jump + rand*jump;
                J(end+1) = j;
            end
        case 'regular'
            % nrepeat*jump = r;
            J = (0:p-1)*jump;
        otherwise
            error('Unknown');
    end
    
    J = mod(round(J),r)+1;
    col = C(J,:);
    
    Q = max(1,min(floor(t*p)+1,p));
    Col = col( Q, :);
    
end

if not(isempty(err))
    %%% Mode error %%%
    err = max(min(err,1),0); % clamp
    q = 512;
    C = jet(q);
    Col = C(floor(err*(q-1))+1,:);
end



hold on;
for i=1:n
    j = mod(i,n)+1;    
    h = plot(gamma([i j]) + eps*(1+1i) );
    set(h,'Color', Col(i,:));
    set(h,'LineWidth', lw);
end
hold off;