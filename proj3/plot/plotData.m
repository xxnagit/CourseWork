function plotData(X1, y1, X2, y2, X3, y3, X4, y4)
%PLOTDATA Plots the data points x and y into a new figure 
%   PLOTDATA(x,y) plots the data points and gives the figure axes labels of
%   population and profit.

% ====================== YOUR CODE HERE ======================
% Instructions: Plot the training data into a figure using the 
%               "figure" and "plot" commands. Set the axes labels using
%               the "xlabel" and "ylabel" commands. Assume the 
%               population and revenue data have been passed in
%               as the x and y arguments of this function.
%
% Hint: You can use the 'rx' option with plot to have the markers
%       appear as red crosses. Furthermore, you can make the
%       markers larger by using plot(..., 'rx', 'MarkerSize', 10);

figure; % open a new figure window
p = plot(X1,y1,X2,y2,X3,y3,X4,y4);
p(1).Marker = '*';
p(2).Marker = '*';
p(3).Marker = '*';
p(4).Marker = '*';

%plot(X, y,'rx','MarkerSize', 5,'-'); %Plot the data
ylabel('Time Spent - s'); %Set the y-axis label
xlabel('Matrix Row/Column Size n'); %Set the x-axis label









% ============================================================

end
