function [Energy,Grad] = compute_rkhs_energy(A,B,options)

% compute_gradient computes the gradient of the enrgy E 

n = length(A);

%%
% Parameters

sigma = getoptions(options, 'sigma', .1);
delta = getoptions(options, 'delta', .1);
rho = getoptions(options, 'rho', .1);

%%
% Kernel (Sum of two gaussian kernels with variance sigma^2 and delta^2)
% Grad is the derivative with respect to the first variable
        
K = @(x,y)exp(-abs(x-y).^2/(2*sigma^2))+ exp(-abs(x-y).^2/(2*delta^2));
GradK = @(x,y)-((x-y).*exp(-abs(x-y).^2/(2*sigma^2)))./sigma^2 -((x-y).*exp(-abs(x-y).^2/(2*delta^2)))./delta^2; 
 

%%
% Define internals callbacks.

%
FwdDer = @(c)(c([2:end 1]) - c);
BwdDer = @(c)(c - c([end 1:end-1]));
CntDer = @(c) FwdDer(c)+ BwdDer(c);

% dot product between two vector fields
dotp = @(a,b)real(a.*conj(b));
Dotp = @(a) sqrt(sum(dotp(a,a)));

% callback for curve derivatives

normalize = @(c)c./max(sqrt(sum(dotp(c,c))),1e-10);
tangent = @(c)normalize(FwdDer(c));
normal = @(c)1i*tangent(c);
curveLength = @(c)sum(abs(FwdDer(c)))/n;

%%
% Tangent, normal and curvature to A.

% lengths
LA = curveLength(A);
LB = curveLength(B);
% tangent
tA = tangent(A);
tB = tangent(B);
% normal
nA = 1i*tA;
nB = 1i*tB;
% normal times curvature
nACurv = BwdDer(tA);
% curvature
CurvA = dotp(nACurv,nA);

%%
% Evaluate the basics quantities.

[J,I] = meshgrid(1:n,1:n);

prevent = @(x)max(x,1e-10);

        
        KAB = K(A(I),B(J));
        KAA = K(A(I),A(J));
        KBB = K(B(I),B(J));
        
        GradKAA = GradK(A(I),A(J));
        GradKAB = GradK(A(I),B(J));
        
        FwdA = FwdDer(A);
        FwdB = FwdDer(B);
        TAB = dotp(FwdA(I),FwdB(J)); % nA.*nB*|A'|*|B'|=tA.*tB
        TAA = dotp(FwdA(I),FwdA(J));
        TBB = dotp(FwdB(I),FwdB(J));
        
        Trap = @(M) 1/4.*(M+M(1:n,[2:n 1])+M([2:n 1],1:n)+M([2:n 1],[2:n 1]));
        
        EAA =  sum(sum(TAA.*Trap(KAA),2),1);
        EBB =  sum(sum(TBB.*Trap(KBB),2),1);
        EAB =  sum(sum(TAB.*Trap(KAB),2),1);
        
        % Distance and Energy
        
        H = EAA + EBB -2*EAB;
        
        Energy = 1/2*H;
      
        
        % C-Gradient
        
        CntA = CntDer(A); 
        CAB = dotp(CntA(I),FwdB(J));
        CAA = dotp(CntA(I),FwdA(J));
        
        GradAB1 = sum(CAB.*(GradKAB + GradKAB(1:n,[2:n 1])),2);
        GradAB2 = ((KAB([n 1:(n-1)],1:n) + KAB([n 1:(n-1)],[2:n 1])- KAB([2:n 1],1:n) - KAB([2:n 1],[2:n 1])))*FwdB;
        GradAB = GradAB1 + GradAB2;
        
        GradAA1 = sum(CAA.*(GradKAA + GradKAA(1:n,[2:n 1])),2);
        GradAA2 = ((KAA([n 1:(n-1)],1:n) + KAA([n 1:(n-1)],[2:n 1])- KAA([2:n 1],1:n) - KAA([2:n 1],[2:n 1])))*FwdA;
        GradAA = GradAA1 + GradAA2;
        

        Grad= (1/4).*(GradAA - GradAB);

     
        % Rigidity matrix
        
        Mcent = diag(2*(abs(FwdA)+ abs(FwdA([n 1:(n-1)]))));
        Mup = diag(abs(FwdA([n 1:(n-2)])),1) + diag(abs(FwdA(n-1)),n-1) ;
        Mdw = diag(abs(FwdA([n 1:(n-2)])),-1) + diag(abs(FwdA(n-1)),-(n-1));
        
        M = (1/(6*n)).*(Mcent +Mup + Mdw);
        
     %   Ai = abs(BwdDer(A));
        
     %   AI = abs(A([3:n 1:2]) - A([2:n 1]));
        
      %  Mcent = diag(1/3*[Ai+AI]+2/3*abs(FwdA));
      %  Mup = diag(1/6*(abs(FwdA([1:(n-1)])) + abs(FwdA([n 1:(n-2)]))),1)  + ...
      %        + diag(1/6*(abs(FwdA(n)) + abs(FwdA(n-1))),n-1) ;
      %  Mdw = diag(1/6*(abs(FwdA([1:(n-1)])) + abs(FwdA([n 1:(n-2)]))),-1) + ...
      %        + diag(1/6*(abs(FwdA(n)) + abs(FwdA(n-1))),-(n-1)) ;
      %  M = 1/n.*(Mcent +Mup + Mdw);
        
        % Gradient with respect to the l^2(A)-inner product
        
        Grad= inv(M)*Grad;            
         

  

end



