%% Machine Learning Class - Exercise 1: Linear Regression

%  Instructions
%  ------------
% 
%  This file contains code that helps you get started on the
%  linear exercise. You will need to complete the following functions 
%  in this exericse:
%
%     plotData.m
%     gradientDescent.m
%     computeCost.m
%
%  For this exercise, you will not need to change any code in this file,
%  or any other files other than those mentioned above.
%
% x refers to the population size in 10,000s
% y refers to the profit in $10,000s
%

%% Initialization
clear ; close all; clc

%% ======================= Part 1: Plotting =======================
fprintf('Plotting Data ...\n')
data = load('data-GFLOP.txt');
X1 = data(1:4,1); y1 = data(1:4,2);
X2 = data(5:8,1); y2 = data(5:8,2);
X3 = data(9:12,1); y3 = data(9:12,2);
X4 = data(13:16,1); y4 = data(13:16,2);
X5 = data(17:20,1); y5 = data(17:20,2);
X6 = data(21:24,1); y6 = data(21:24,2);
X7 = data(25:28,1); y7 = data(25:28,2);


%m = length(y); % number of training examples

% Plot Data
% Note: You have to complete the code in plotData.m
plotData(X1, y1, X2, y2, X3, y3, X4, y4);
plotData(X5, y5, X6, y6, X7, y7, X7, y7);

fprintf('Program paused. Press enter to continue.\n');
pause;

