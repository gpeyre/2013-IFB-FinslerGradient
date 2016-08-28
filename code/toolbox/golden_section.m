function tau = golden_section( func, a, b, niter )

alpha= (sqrt(5)-1)/2;

for i=1:niter
    
    c = alpha*a + (1-alpha)*b;
    d= a+b-c;
    
    if func(c)<func(d)
        b=d;
    else
        a=c;
    end    

tau= c;

end

