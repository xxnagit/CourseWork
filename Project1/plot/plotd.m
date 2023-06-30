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
data = load('data-mac.txt');
X1 = data(1:12, 1); y1 = data(1:12, 2);
X2 = data(13:24,1); y2 = data(13:24,2);
X3 = data(25:36,1); y3 = data(25:36,2);
X4 = data(37:48,1); y4 = data(37:48,2);
X5 = data(49:60,1); y5 = data(49:60,2);
X6 = data(61:72,1); y6 = data(61:72,2);
X7 = data(73:84,1); y7 = data(73:84,2);
X8 = data(85:96,1); y8 = data(85:96,2);

%m = length(y); % number of training examples

% Plot Data
% Note: You have to complete the code in plotData.m
plotData(X1, y1, X2, y2, X3, y3);
plotData(X3, y3, X4, y4, X5, y5);
plotData(X3, y3, X4, y4, X5, y5);

fprintf('Program paused. Press enter to continue.\n');
pause;

