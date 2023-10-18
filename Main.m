% This is the matlab code for the optimization algorithm, namely Bonobo Optimizer (BO).
% This is written for solving unconstrained optimization problems. However, it can also solve constrained optimization
% problems with penalty function approaches. 
% Moreover, this for solving minimization problems.
% For details of the BO algorithm, kindly refer and cite as mentioned below:
% A. K. Das and D. K. Pratihar, "Bonobo optimizer (BO): an intelligent heuristic with selfadjusting parameters over continuous spaces and its applications to engineering problems," 
% Applied Intelligence, 2021, DOI: 10.1007/s10489-021-02444-w
% For any query, please email to: amit.besus@gmail.com
clc;close all;clear all;
tic;   % CPU time measure
CostFunction = @(x)MyObjectiveFunction(x); % Objective function 
d=4;  % No. of Variables
Var_min=[-100 -100 -100 -100];  % Lower variable Boundaries
Var_max=[100 100 100 100];   % Upper variable Boundaries
%% Common parameters of  BO similar to  other optimization algorithms
N=30; % No. of bonobos in the population, i.e. population size
max_it=100;  % Maximum number of iterations
[bestcost,alphabonobo,convergence_curve]=BO(N,d,Var_min,Var_max,max_it,CostFunction);
disp(['Best Cost: ' num2str(bestcost)]);
disp(['Bestsolution: ' num2str(alphabonobo)]);
figure
plot (1:max_it,convergence_curve,'-*')
title('Convergence Curve')
xlabel('Number of iterations')
ylabel('Evolution of best objective value')
toc;