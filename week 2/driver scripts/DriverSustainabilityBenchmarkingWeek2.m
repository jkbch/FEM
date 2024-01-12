% DriverSustainabilityBenchmarkingWeek2.m
%
% This is a driver script that will be used for assessment of
% week 2 handins in DTU course 02623 Finite Element Method for Partial
% Differential Equations.
%
% Participants of the course supplies the following drivers.
% Test that they work by running this code before submission)
%   - Driver28b.m
%   - Driver28c.m
%
% Developed by Allan P. Engsig-Karup.
close all
clear all
clc

%% TODO PUT IN YOUR GRROUP NO AND STUDENT IDs FOR THE GROUP HERE
groupNo = <X>; % Put your group no here.
groupStudentIDs = '<Y/Z>'; % Put your student id's here.

%% PATHS 
cd(fileparts(which(mfilename))) % change to directory
dirThisScript   = cd;           % store path
dirStoreResults = '/Users/apek/02623/Handinresults/';

%% PARAMETERS FOR SUUSTAINABILITY CALCULATION
CO2intensity  = 0.285; % [kg CO2/kWh], https://communitiesforfuture.org/collaborate/electricity-map/
PowerEstimate = 60; % [kW]

%% PARAMETERS FOR THE SELECTED EXERCISES (DO NOT CHANGE)
%% Define input parameters
x0 = 0; y0 = 0;
L1 = 1; L2 = 1;
noelms1=40;
noelms2=50;
lam1 = 1;
lam2 = 1;
fun = @(x,y) cos(pi*x).*cos(pi*y);
qt  = @(x,y) 2*pi^2*cos(pi*x).*cos(pi*y);

%% EXECUTE CODE
% Let's call the FEM BVP 2D Solver you produced
% time the code  using tic and toc

%% Call Group 30 solver
tic;
[VX,VY,EToV,U] = Driver28b(x0,y0,L1,L2,noelms1,noelms2,lam1,lam2,fun,qt);
tend = toc;
DOF1 = length(U(:));


x0 = -1; y0 = -1;
L1 = 2; L2 = 2;
tic;
[VX2,VY2,EToV2,U2] = Driver28c(x0,y0,L1,L2,noelms1,noelms2,lam1,lam2,fun,qt);
tend2 = toc;
DOF2 = length(U2(:));

CPUtime1 = tend;
CPUtime2 = tend2;
CO2eq1 = CPUtime1/3600*PowerEstimate/1000*CO2intensity; 
CO2eq2 = CPUtime2/3600*PowerEstimate/1000*CO2intensity; 

%% Visualization
figure
subplot(1,2,1)
trisurf(EToV,VX,VY,U)
xlabel('x')
ylabel('y')
zlabel('u(x,y)')
title(sprintf('2.8b. Group: %s, Time: %.4e, DOF: %d, noelsm1=%d, noelms2=%d, CO2e=%.4e',groupStudentIDs,tend,DOF1,noelms1,noelms2,CO2eq1))

%figure
subplot(1,2,2)
trisurf(EToV2,VX2,VY2,U2)
xlabel('x')
ylabel('y')
zlabel('u(x,y)')
title(sprintf('2.8c. Group: %s, Time: %.4e, DOF: %d, noelsm1=%d, noelms2=%d, CO2e=%.4e',groupStudentIDs,tend2,DOF2,noelms1,noelms2,CO2eq2))

set(gcf,'position',[127         372        1307         576])

%% STORE THE RESULTS
filename=sprintf('Week2ResultsGroup%d.txt',groupNo); 
fid=fopen([dirStoreResults filename],'w'); % open file identifier
formatSpec = '%d \n';
fprintf(fid, formatSpec, groupNo);
formatSpec = '%s \n';
fprintf(fid, formatSpec, groupStudentIDs);
formatSpec = '%.4e %d .4e \n';
fprintf(fid, formatSpec, [CPUtime1, DOF1, CO2eq1]);
fprintf(fid, formatSpec, [CPUtime2, DOF2, CO2eq2]);
fclose(fid) %close file indentifier


