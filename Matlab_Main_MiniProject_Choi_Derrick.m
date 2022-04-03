%% ASEN 5007 - Mini Project - Main
% Purpose: develop a Bubnov-Galerkin finite element approximation to
% variational form of 1D-model problem and use it solve a problem on heat
% conduction and a cylindrical rod undergoing a tensile test
%
% Author: Derrick Choi
% Date Created: 9/26/2021
% Last Modified: 10/14/2021
clear;clc;close all;
%% Verification Test 1 Homework Problem
fprintf('Problem 2\n')
fprintf('Performing Verification Test 1: with Homework 1 Problem 1c and 1d\n\n')
%Test 1: Check if Homework Problem 1 is recovered for k = 3 and k = 1 Bubnov-Galerkin
kappa = @(x) -1;
f = @(x) x;
g_0 = 0;
g_L = 0;
L = 1;
%exact solutions
u_HW1c = @(x) x/6.*(x.^2-1);
uprime_HW1c = @(x) 1/6*(3.*x.^2-1);
[x_hw1c,d_hw1c] = One_Dim_Model_Problem(1,3,kappa,f,g_0,g_L,L);
[x_hw1d,d_hw1d] = One_Dim_Model_Problem(3,1,kappa,f,g_0,g_L,L);
H1_error_1d = H1NormError(u_HW1c,uprime_HW1c,d_hw1d,3,1,L,g_0,g_L);
fprintf('H1-norm error of HW 1 Problem 1d = %0.3e\n',H1_error_1d); 

% Finite element approximations
n_el = ceil(logspace(1,3));
errorsk1 = zeros(1,length(n_el));
errorsk2 = zeros(1,length(n_el));
errorsk3 = zeros(1,length(n_el));
for e = 1:length(n_el)
    [~,d_k1] = One_Dim_Model_Problem(1,n_el(e),kappa,f,g_0,g_L,L);
    [~,d_k2] = One_Dim_Model_Problem(2,n_el(e),kappa,f,g_0,g_L,L);    
    [~,d_k3] = One_Dim_Model_Problem(3,n_el(e),kappa,f,g_0,g_L,L);
    errorsk1(e) = H1NormError(u_HW1c,uprime_HW1c,d_k1,1,n_el(e),L,g_0,g_L);
    errorsk2(e) = H1NormError(u_HW1c,uprime_HW1c,d_k2,2,n_el(e),L,g_0,g_L);
    errorsk3(e) = H1NormError(u_HW1c,uprime_HW1c,d_k3,3,n_el(e),L,g_0,g_L);
end
fprintf('Plotting H1-norm error vs number of elements for Verification Test 1\n\n')
figure('Name','Rates of Convergence Test 1: HW1 Problem 1c and 1d','Position',[400 100 800 600])
loglog(n_el,errorsk1,'LineWidth',1.5)
hold on
loglog(n_el,errorsk2,'LineWidth',1.5)
loglog(n_el,errorsk3,'LineWidth',1.5)
set(gca,'TickLabelInterpreter','latex','FontSize',12)
xlabel('Number of elements','Interpreter','latex','FontSize',14)
ylabel('$\|u-u^h\|_{H^1}$','Interpreter','latex','FontSize',14)
title('$H^1$ Norm Error vs Number of Elements ($u(x) = x/6(x^2-1)$, $\kappa(x) = 1$, $x \in (0,1)$)','FontSize',16,'Interpreter','latex')
legend('k = 1','k = 2','k = 3','Interpreter','latex','FontSize',12)
grid on
%% Verification Test 2 Maufactured Solution

fprintf('Performing Verification Test 2: Manufactured solution u = exp(-x)sin(pi*x)\n')
%problem set up
kappa = @(x) 1;
f = @(x) (pi^2-1)*exp(-x)*sin(pi*x)+2*pi*exp(-x)*cos(pi*x);
L = 1;
u_M1 = @(x) exp(-x)*sin(pi*x);
u_prime_M1 = @(x) exp(-x)*(pi*cos(pi*x)-sin(pi*x));
g0 = 0;
gL = 0;

errorsk1 = zeros(1,length(n_el));
errorsk2 = zeros(1,length(n_el));
errorsk3 = zeros(1,length(n_el));

for e = 1:length(n_el)
    [~,d_k1] = One_Dim_Model_Problem(1,n_el(e),kappa,f,g0,gL,L);
    [~,d_k2] = One_Dim_Model_Problem(2,n_el(e),kappa,f,g0,gL,L);
    [~,d_k3] = One_Dim_Model_Problem(3,n_el(e),kappa,f,g0,gL,L);
    errorsk1(e) = H1NormError(u_M1,u_prime_M1,d_k1,1,n_el(e),L,g0,gL);
    errorsk2(e) = H1NormError(u_M1,u_prime_M1,d_k2,2,n_el(e),L,g0,gL);
    errorsk3(e) = H1NormError(u_M1,u_prime_M1,d_k3,3,n_el(e),L,g0,gL);    
end

fprintf('Plotting H1-norm error vs number of elements for Verification Test 2\n')
figure('Name','Rates of Convergence Test 2: u = exp(-x)sin(x)','Position',[400 100 800 600])
loglog(n_el,errorsk1,'LineWidth',1.5)
hold on
loglog(n_el,errorsk2,'LineWidth',1.5)
loglog(n_el,errorsk3,'LineWidth',1.5)
set(gca,'TickLabelInterpreter','latex','FontSize',12)
xlabel('Number of elements','Interpreter','latex','FontSize',14)
ylabel('$\|u-u^h\|_{H^1}$','Interpreter','latex','FontSize',14)
title('$H^1$ Norm Error vs Number of Elements ($u = e^{-x}sin(\pi x)$, $\kappa(x) = 1$, $x\in (0,1)$)','FontSize',16,'Interpreter','latex')
legend('k = 1','k = 2','k = 3','Interpreter','latex','FontSize',12)
grid on
%% Verification Test 3 Manufacture Solution 2
fprintf('\nPerforming Verification Test 3: Manufactured Solution u(x) = exp(x), kappa(x) = x\n')
%problem set-up
kappa = @(x) x;
u_M2 = @(x) exp(x);
u_prime_M2 = @(x) exp(x);
L = 1;
g0 = 1;
gL = exp(1);
f = @(x) -exp(x)*(1+x);

errorsk1 = zeros(1,length(n_el));
errorsk2 = zeros(1,length(n_el));
errorsk3 = zeros(1,length(n_el));

for e = 1:length(n_el)
    [~,d_k1] = One_Dim_Model_Problem(1,n_el(e),kappa,f,g0,gL,L);
    [~,d_k2] = One_Dim_Model_Problem(2,n_el(e),kappa,f,g0,gL,L);
    [~,d_k3] = One_Dim_Model_Problem(3,n_el(e),kappa,f,g0,gL,L);
    errorsk1(e) = H1NormError(u_M2,u_prime_M2,d_k1,1,n_el(e),L,g0,gL);
    errorsk2(e) = H1NormError(u_M2,u_prime_M2,d_k2,2,n_el(e),L,g0,gL);
    errorsk3(e) = H1NormError(u_M2,u_prime_M2,d_k3,3,n_el(e),L,g0,gL);    
end

fprintf('Plotting H1-norm error vs number of elements for Verification Test 3\n')
figure('Name','Rates of Convergence Test 3: u = exp(x)','Position',[400 100 800 600])
loglog(n_el,errorsk1,'LineWidth',1.5)
hold on
loglog(n_el,errorsk2,'LineWidth',1.5)
loglog(n_el,errorsk3,'LineWidth',1.5)
set(gca,'TickLabelInterpreter','latex','FontSize',12)
xlabel('Number of elements','Interpreter','latex','FontSize',14)
ylabel('$\|u-u^h\|_{H^1}$','Interpreter','latex','FontSize',14)
title('$H^1$ Norm Error vs Number of Elements ($u = e^{x}$, $\kappa(x) = x$, $x\in (0,1)$)','FontSize',16,'Interpreter','latex')
legend('k = 1','k = 2','k = 3','Interpreter','latex','FontSize',12)
grid on

%% Problem 3 - Heat problem
fprintf('\nProblem 3\n')
%Problem set-up
R = 0.025;
L = 2;
Tend = 30;
h_max = 1500;
kappa = @(x) 385;
h = @(x) h_max*exp(-100*(x/L-0.5)^2);
f = @(x) 2/R*h(x);

%solve problem with one element
[x_T_n1,d_T_n1] = One_Dim_Model_Problem(1,2,kappa,f,Tend,Tend,L);
maxT_n1 = max([Tend;d_T_n1;Tend]);

%solve problem with two elements
[x_T_n2,d_T_n2] = One_Dim_Model_Problem(3,10,kappa,f,Tend,Tend,L);
maxT_n2 = max([Tend;d_T_n2;Tend]);

%keep iterating with more elements until relative change between max
%temperatures is small
n = 20;
while abs((maxT_n2-maxT_n1)/maxT_n2) > 1e-6

    %set old solution to current solution
    x_T_n1 = x_T_n2;
    d_T_n1 = d_T_n2;
    maxT_n1 = maxT_n2;
    
    %find new solution with more elements
    [x_T_n2,d_T_n2] = One_Dim_Model_Problem(3,n,kappa,f,Tend,Tend,L);
    maxT_n2 = max([Tend;d_T_n2;Tend]);
    %update number of elements used
    n = n+10;
end
fprintf('Max Temperature = %f degrees C\n',maxT_n2)
%% Problem 4 - Tensile Test
fprintf('\nProblem 4\n')
E = 200e9;
Amax = 1e-4;
Amin = 2e-5;
u0 = 0;
u_end = 2e-5;
L = 0.1;

A = @(x) Amax-(Amax-Amin)*exp(-50*(x/L-0.5)^2);

kappa = @(x) E*A(x);
f = @(x) 0;

%solve problem with one elements
[u1,d_u1] = One_Dim_Model_Problem(3,1,kappa,f,u0,u_end,L);
dudx1 = diff([u0;d_u1;u_end])./diff(u1');
max_strain1 = max(dudx1);

%solve problem with two elements
[u2,d_u2] = One_Dim_Model_Problem(3,100,kappa,f,u0,u_end,L);
dudx2 = diff([u0;d_u2;u_end])./diff(u2');
max_strain2 = max(dudx2);

%keep solving with more elements until relative change in max derivative is small
n = 200;
while abs((max_strain2-max_strain1)/max_strain2) >1e-6
    
    %set the old solution to current solution
    u1 = u2;
    d_u1 = d_u2;
    dudx1 = dudx2;
    max_strain1 = max_strain2;
    
    %find new solution with more elements
    [u2,d_u2] = One_Dim_Model_Problem(3,n,kappa,f,u0,u_end,L);
    dudx2 = diff([u0;d_u2;u_end])./diff(u2');
    max_strain2 = max(dudx2);
    
    %update number of elements
    n = n+100;
end
%find max stress and output to command window
max_axial_stress = E*max_strain2;
fprintf('Max Axial Stress = %0.2f Pa \n',max_axial_stress)
