
clear all;
close all;
clc;

%% Specify file for data extraction
open('fit_byom_metabolism_clomazone_1.fig') ; % choose corresponding file name of MATLAB figure
h = findobj(gcf, 'Type', 'line');

% NOTE: As there were no input data for the metabolite provided at 07
% degrees, the locations of the data for the plotted lines within the
% matlab figure are different. Make sure to use the correct locations by
% identifying them first (i.e., open xdata and ydata, first lines in
% each section)

%% Extract xdata
xdata= get(h,'XData'); % open variables stored in XData for identification
xdata1= get(h(11,1),'XData'); % X data for all lines

% %% Find out Handles:
%  numHandles = length(h);
% disp(numHandles);  % Display the number of lines found
%     11
% 
% for i = 1:numHandles
%     fprintf('Handle %d:\n', i);
%     disp(get(h(i), 'YData'));  % Display YData for each line
%end
%% Extract ydata
ydata= get(h,'YData'); % open variables stored in YData for identification
% % Total internal conentrations (in µmol/kg) 
% NOTE: Nominal exposure concentrations are marking the lines
ydata1= get(h(9,1),'YData'); % mean              
ydata2= get(h(7,1),'YData'); % upper Limit       
ydata3= get(h(8,1),'YData'); % lower limit       

% % Internal conentrations in receptor-complex/membrane (in µmol/kg) 
% NOTE: Nominal exposure concentrations are marking the lines
ydata4= get(h(6,1),'YData'); % mean              
ydata5= get(h(4,1),'YData'); % upper Limit       
ydata6= get(h(5,1),'YData'); % lower limit       

% % Internal conentrations supertenant+debris (in µmol/kg) 
% NOTE: Nominal exposure concentrations are marking the lines
ydata7= get(h(12,1),'YData'); % mean             
ydata8= get(h(10,1),'YData'); % upper Limit       
ydata9= get(h(11,1),'YData'); % lower limit      

%% Export data to txt
fig01 = []; %create empty table
fig01(:,1) = xdata1 ; % xdata for all lines
fig01(:,2) = ydata1 ; 
fig01(:,3) = ydata2 ;
fig01(:,4) = ydata3 ;
fig01(:,5) = ydata4 ; 
fig01(:,6) = ydata5 ;
fig01(:,7) = ydata6 ;
fig01(:,8) = ydata7 ; 
fig01(:,9) = ydata8 ;
fig01(:,10) = ydata9 ;

dlmwrite('fit_byom_metabolism_clomazone_1.txt', fig01, ','); % write dataframe in txt
