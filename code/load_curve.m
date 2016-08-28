function gamma = load_curve(name, n)

% This function loads the several curves used in the numerical tests

switch name

    case 'rod'
        gamma = [0 1 1 0]' + 1i*[0 0 1/10 1/10]';
        
    case 'harm1'
        gamma = [0 1/2 1/2 1 1 1/2 1/2 0]' + 1i*[0 0 3/5 2/5 3/5 4/5 1 1]';
    
    case 'harm2'
        gamma = [0 1/2 1/2 3/4 1 1 3/4  1/2 1/2 0]' + 1i*[0 0 3/5 1/2 3/5 4/5 7/10 4/5 1 1]';       
 
    case  {'horse1' 'horse2'}
        f = load_image(name);
        gamma = perform_levelset_extraction(f,n);     
                     
    case  {'man1' 'man2'}
        f = load_image(name);
        gamma = perform_levelset_extraction(f,n);  
        
                  
end

isper =1;
gamma = perform_curve_interpolation(gamma,n, [], isper);
