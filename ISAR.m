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
r = 170 * millimeters;

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
for m = 1 : M
    for n = 1 : N
        t(m,n) = (n-1)*T2 + (m-1)*(N-1)*T2;
    end
end

% DETERMINE SIG IN LINEAR SCALE
SIG = 10^(SIGdBsm/10) * (meters)^2;

% INITIALIZE TRANSFER FUNCTION ARRAYS
HF = zeros(M,N);

%
% MAIN LOOP
%

% POPULATE TRANSFER MATRIX
for m = 1 : M
    for n = 1 : N
        HF(m,n) = sqrt(SIG) * (exp(1i*4*pi*f0(n)*R0/c0) + ...
        exp((1i*4*pi*f0(n)*(R0 + r*cos((pi/4) + OMEGA*t(m,n))))/c0) + ...
        exp((1i*4*pi*f0(n)*(R0 + r*cos((3*pi/4) + OMEGA*t(m,n))))/c0) + ...
        exp((1i*4*pi*f0(n)*(R0 + r*cos((5*pi/4) + OMEGA*t(m,n))))/c0) + ...
        exp((1i*4*pi*f0(n)*(R0 + r*cos((7*pi/4) + OMEGA*t(m,n)))/c0)));
    end
end

% PERFORM 2D FFTs
HF = fftshift(fft(HF,N,1)/N);   % Column-Wise FFT
HF = fft(HF,M,2)/M;             % Row-Wise FFT

% CONVERT TO dBsm
HF = 20*log10(abs(HF));         % 20 used because it is a voltage ratio

% PLOT RESULTING IMAGE
xa = [0 : N-1]*dRd; xa = xa - mean(xa); % Generate x axis
ya = [0 : M-1]*dRc; ya = ya - mean(ya); % Generate y axis
h = imagesc(xa,ya,HF);
h = get(h,'Parent');
set(h,'FontSize',11,'YDir','normal');
axis equal tight
colormap('Jet');
c = colorbar;
caxis([-60 -10]);
title('$\textrm{ISAR Image of a Drone}$','Interpreter','LaTex',...
    'FontSize',14);
xlabel('$\textrm{x (m)}$','Interpreter','LaTex','FontSize',12);
ylabel('$\textrm{y (m)}$','Interpreter','LaTex','FontSize',12,...
    'Rotation',90);

