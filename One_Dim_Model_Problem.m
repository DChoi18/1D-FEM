function [x,d] = One_Dim_Model_Problem(k,n_el,kappa,f,g_0,g_L,L)
% function that approximately solves the one dimensional model problem
% -d/dx(k(x)du/dx) = f(x) using a Bubnov-Galerkin finite element
% approximation of the variational form.
%
% Need to fix: use the derivative of the basis function
%
% Inputs:
%         k = polynomial degree
%         n_el = number of elements
%         kappa = function handle (material modulus)
%         f = function handle (forcing term)
%         g_0 = boundary condition at x = 0
%         g_L = boundary condition at x = L
%         L = length of domain
% Outputs:
%         x = node locations
%         d = finite element solutions at the nodes
% Author: Derrick Choi
%% Set up

% Determine number of basis functions
n = (k-1)*n_el+n_el-1;
K = zeros(n,n);
F = zeros(n,1);

h_e = L/n_el; %element size
x = 0:h_e:L; %element node locations
%% Element Assembly
for e = 1:n_el %loop over elements
    % Get element stiffness and force
    [k_e,f_e] = Element_K_F(e,k,kappa,f,h_e,L);
    for a = 1:k+1
        A = k*(e-1)+a-1;
        if A<=n && A >=1
            for b = 1:k+1
                B = k*(e-1)+b-1;
                if B<=n && B >=1
                    K(A,B) = K(A,B)+k_e(a,b);
                elseif B ==0
                    F(A) = F(A)-k_e(a,b)*g_0;
                elseif B == n+1
                    F(A) = F(A)-k_e(a,b)*g_L;
                end
            end
            F(A) = F(A)+f_e(a);
        end
    end
end
%% Solve Kd = F
d = K\F;
%% Output x for corresponding DOFs and the endpoints
x = 0:L/(n_el*(k-1)+n_el):L;
end

function [k_e,f_e] = Element_K_F(A,k,kappa,f,h_e,L)
%% Element Formation
k_e = zeros(k+1,k+1);
f_e = zeros(k+1,1);

x = 0:h_e:L; %node locations

%determine quadrature points and weights for Gaussian quadrature 
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

xi = -1:2/k:1; %node locations in parent element
for q = 1:length(xi_q)
    for a = 1:k+1
        xe_a = x(A); %local node location
        xe = @(y) xe_a+h_e*((y+1)/2); %inverse mapping function back to physical element
        
        point = xe(xi_q(q)); %quadrature point location on physical element
        kappa_point = kappa(point); %kappa value at quadrature point
        f_point = f(point);         %f value at quadrature point
        
        %define parent basis function and its derivative for 'a' loop
        Nhat_a = @(xi_point) prod((xi(xi~=xi(a))-xi_point)./(xi(xi~=xi(a))-xi(a)));
        w = prod(1./(xi(xi~=xi(a))-xi(a)));
        s = @(x) 0;
        for i = 1:k
            xi_not_a = xi(xi~=xi(a));
            s = @(x) s(x) + prod( xi( xi~=xi(a) & xi ~= xi_not_a(i) )-x );
        end
        dNhat_a = @(xi_point) w*s(xi_point);
        for b = 1:k+1
            %define derivative of parent basis function for 'b' loop
            wb = prod(1./(xi(xi~=xi(b))-xi(b)));
            s2 = @(x) 0;
            for j = 1:k
                xi_not_b = xi(xi~=xi(b));
                s2 = @(x) s2(x) + prod(xi(xi~=xi(b) & xi ~= xi_not_b(j))-x); 
            end
            dNhat_b = @(x) wb*s2(x);
            %add contribution to integral at quadrature point
            k_e(a,b) = k_e(a,b)+2/h_e*kappa_point*dNhat_a(xi_q(q))*dNhat_b(xi_q(q))*wq(q);
        end
        f_e(a) = f_e(a)+h_e/2*f_point*Nhat_a(xi_q(q))*wq(q);
    end
end

end