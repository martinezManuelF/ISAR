% =========================================================================
% FILE DESCRIPTION
% =========================================================================
% ISAR.m
% Inverse Synthetic Aperture Radar Simulation
%
% =========================================================================
% TEAM MEMBERS
% =========================================================================
% Herrera, Cesar
% Martinez, Manuel
% Ontiveros, Raul
% Salais, Irvin
%
% =========================================================================
% COURSE
% =========================================================================
% EE5389-Radar Signal Processing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RESTORE STATE
clear all;
close all;
clc;

% UNITS
meters       = 1;
centimeters  = 1e-2 * meters;
millimeters  = 1e-3 * meters;
seconds      = 1;
milliseconds = 1e-3 * meters;
microseconds = 1e-6 * meters;
hertz        = 1/seconds;
kilohertz    = 1e3 * hertz;
megahertz    = 1e6 * hertz;
gigahertz    = 1e9 * hertz;
degrees      = pi/180;

% CONSTANTS
c0 = 3e8 * meters/seconds;

% OPEN FIGURE WINDOW
figure('Color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RADAR CONSTRAINTS
f1  = 30 * gigahertz;    % Initial frequency
T2  = 1 * microseconds;  % PRT
R0  = 100 * meters;      % Distance from center point target to radar
dRd = 20 * millimeters;  % Down range resolution
dRc = 20 * millimeters;  % Cross range resolution

% TARGET CONSTRAINTS
SIGdBsm = -10;
r = 30 * millimeters;

% IMAGE CONSTRAINTS
Sx = 640 * millimeters;
Sy = 640 * millimeters;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DETERMINE RADAR BANDWIDTH
BETA = c0/(2*dRd);

% DETERMINE NUMBER OF FREQUENCY STEPS (N)
N = Sx/dRd;

% DETERMINE NUMBER OF PULSE TRAINS (M)
M = Sy/dRc;

% DETERMINE FREQUENCY STEP
df = BETA/(N-1);

% DETERMINE INTEGRATION TIME
Tint = M*(N-1)*T2;

% DETERMINE FINAL FREQUENCY
fmax = f1 + N*df;

% DETERMINE CENTER FREQUENCY AND CENTER WAVELENGTH
fcenter = (f1 + fmax)/2;
lamc    = c0/fcenter;

% DETERMINE RANTE OF CHANGE OF ANGLE
dTHETA = lamc/(2*dRc);

% DETERMINE OMEGA
OMEGA = dTHETA/Tint;

% CREATE A FREQUENCY ARRAY
f0 = linspace(f1,fmax,N);

% CREATE THE TIME ARRAY
for n = 1 : N
    for m = 1 : M
        t(m,n) = (n-1)*T2 + (m-1)*(N-1)*T2
    end
end

%CHANGES
