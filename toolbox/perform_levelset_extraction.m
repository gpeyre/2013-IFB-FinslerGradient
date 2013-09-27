function A = perform_levelset_extraction(f,n,t)

if nargin<3
    f = rescale(sum(f,3));
    t = .5;
end

p = max(size(f,1), size(f,2));
c = contourc(f,[1 1]*t);

lmax = 0;
A = [];
while not(isempty(c))
    lc = c(2,1); c(:,1) = [];
    B = c(:,1:lc);  c(:,1:lc) = [];
    if lc>length(A)
        A = B;
    end
end
B = (B-1)/(p-1);
A = B(2,:) + 1i*B(1,:);
A = A(:);
% A = perform_curve_interpolation(A,n);