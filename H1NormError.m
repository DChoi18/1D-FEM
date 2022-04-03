function error = H1NormError(u,uprime,d,k,nel,L,g0,gL)
% function to compute the H1-error between the analytical solutiona and the
% finite element solution to the 1D-model problem
%
% Inputs:
%         u = analytical solution
%         uprime = derivative of analytical solution
%         d = degrees of freedom solved from Bubnov-Galerkin discretization
%         k = polynomial degree used in finite element solution
%         nel = number of elements
%         L = length of domain
%         g0, gL = boundary conditions at 0 and L respectively
% Outputs: 
%         error  = H1-norm error
% Author: Derrick Choi
h_e = L/nel;
nodes = 0:h_e:L;
xi = -1:2/k:1;
n = (k-1)*nel+nel-1;
%determine quadrature points and weights for Gaussian quadrature (use nq = k+1)
if k == 1
    xi_q = [-1/sqrt(3),1/sqrt(3)];
    wq = [1,1];
elseif k == 2
    xi_q = [0,-sqrt(3/5),sqrt(3/5)];
    wq = [8/9,5/9,5/9];
else
    xi_q = [sqrt(3/7-2/7*sqrt(6/5)),-sqrt(3/7-2/7*sqrt(6/5)),...
        sqrt(3/7+2/7*sqrt(6/5)),-sqrt(3/7+2/7*sqrt(6/5))];
    wq = [(18+sqrt(30))/36,(18+sqrt(30))/36,...
        (18-sqrt(30))/36,(18-sqrt(30))/36];
end

error1 = 0; %first term in H1-norm
error2 = 0; %derivative term in H1-Norm
%compute squared L2 Norms
for e = 1:nel
    % initialize discrete solutions to zero
    uh = @(x) 0;
    uh_prime = @(x) 0;
    for a = 1:k+1
        A = k*(e-1)+a-1; %global node
        %Basis function
        Nhat_a = @(x) prod( (xi(xi~=xi(a))-x )./( xi(xi~=xi(a))-xi(a) ));
        %Form Derivative of basis function
        w = prod(1./(xi(xi~=xi(a))-xi(a)));
        s = @(x) 0;
        for i = 1:k
           xi_not_a =  xi(xi~=xi(a)); % xi values that are excluded from product in the basis function
           s = @(x) s(x)+ prod( xi( xi~=xi(a) & xi ~= xi_not_a(i))-x);
        end
        dNhat_a = @(x) -w*s(x);
        % Add contributions of basis function to uh on the current element
        if A == 0
            uh = @(x) uh(x)+g0*Nhat_a(x);
            uh_prime = @(x) uh_prime(x)+g0*dNhat_a(x);
        elseif A <=n && A >=1
            uh = @(x) uh(x)+d(A)*Nhat_a(x);
            uh_prime = @(x) uh_prime(x)+d(A)*dNhat_a(x);
        elseif A == n+1
            uh = @(x) uh(x)+gL*Nhat_a(x);
            uh_prime = @(x) uh_prime(x) + gL*dNhat_a(x);
        end
    end
    xe_a = nodes(e); %node location in physical space
    xe = @(y) xe_a+h_e*((y+1)/2); %inverse mapping function back to physical element
    
    %sum over quadrature points
    for q = 1:length(xi_q)
        error1 = error1+(u(xe(xi_q(q)))-uh(xi_q(q)))^2*h_e/2*wq(q);
        error2 = error2+(uprime(xe(xi_q(q)))-2/h_e*uh_prime(xi_q(q)))^2*h_e/2*wq(q);
    end
end

%Output the H1 norm error
error = sqrt(error1+error2);
end