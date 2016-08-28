function isok = load_mosek(custum_path)

% load_mosek - add mosek software to the path.
%
%   isok = load_mosek(custum_path);
%
%   Type:
%       mosekopt;
%   to check installation.
%
%   Copyright (c) 2012 Gabriel Peyre

mosek_location = {'/Users/nardi/Dropbox/banach-optim/NumericalTests/mosek/', ...
                  '/Users/gabrielpeyre/Dropbox/matlab/mosek/','/Users/vialard/Dropbox/banach-optim/matlab/mosek/'};
  

if nargin>0
    mosek_location{end+1} = custum_path;
end

isok = 0;
for i=1:length(mosek_location)
    if exist(mosek_location{i})==7
        addpath([mosek_location{i} '6/toolbox/r2009b']);
        setenv('MOSEKLM_LICENSE_FILE', [mosek_location{i}  '/6/licenses/mosek.lic']);
        isok = 1;
        return;
    end    
end

warning('Failed to located Mosek');