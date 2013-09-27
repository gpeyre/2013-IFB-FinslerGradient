function [gamma1,s1] = perform_curve_interpolation(gamma,p,s, isper)

% perform_curve_interpolation - reinterpolate a curve evenly
%
%   [c1,s1] = perform_curve_interpolation(c0,p,[s], [isper]);
%
%   p is the new number of points.
%   s is a bijective mapping from the curve to [0,1], that is updated
%       according to the new parameterization.
%   set isper=1 for periodic (closed) and isper=0 for open curves
%   (if not set, tries to automatically detect it).
%
%   Copyright (c) 2010 Gabriel Peyre

method = 'spline';
method = 'linear';

if nargin<4
    isper=isperiodic(gamma);
end

s1 = [];
if isper
    %%% Periodic curves %%%
    gamma1 = gamma;
    gamma1(end+1) = gamma1(1);
    % new absicse of the points
    d = abs(gamma1(1:end-1)-gamma1(2:end));
    d = max(d,1e-10);
    d = [0;cumsum(d)];
    d = d/d(end);
    % interpolation
    gamma1 = interp1(d,gamma1,(0:p-1)'/p, method); 
    % re-sample the parameterization
    if nargin>2  && not(isempty(s))
        s(end+1) = 1;
        s1 = interp1(d,s,(0:p-1)'/p, method);
    end    
else
    %%% Non-Periodic curves %%%
    % new absicse of the points
   	d = [0;cumsum( 1e-5 + abs(gamma(1:end-1)-gamma(2:end)) )];
    d = d/d(end);
    % interpolation
    gamma1 = interp1(d,gamma,(0:p-1)'/(p-1), method); 
    % re-sample the parameterization
    if nargin>2  && not(isempty(s))
        s1 = interp1(d,s,(0:p-1)'/(p-1), method);
    end
end
