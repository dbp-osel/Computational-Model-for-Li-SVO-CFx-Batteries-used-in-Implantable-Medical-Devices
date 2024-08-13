%% This is a Thermal-Electrochemical simulation model for the Li/SVO-CFx primary battery 
% There exists three state variables (x1= SOC, x2= transient voltage, x3=
% battery temperature). The input to this model is discharge current, and
% the output is OCV.
%The code solves for the states and output with a forward Euler method.
% Code by Mahsa Doosthosseini 2022-2023

clear all
close all
clc

BodyTemperature = 37; % Degree [C] 
ambientT = BodyTemperature+273.15; % [K] ± 10 degree C
load('OCV_DOD_Li_CFxSVO_fine.mat') % Extracted from OCV-DOC curve of Li/SVO-CFx battery from Gomadam et.al 
% [Gomadam, Parthasarathy M., Don R. Merritt, Erik R. Scott, Craig L. Schmidt, Paul M. Skarstad, and John W. Weidner. 
% "Modeling lithium/hybrid-cathode batteries." Journal of power sources 174, no. 2 (2007): 872-876]
ocv = OCV_LiCFxSVO_fine; %Li/CFx-SVO
soc = (100-DOD_LiCFxSVO_fine)/100;
load('th.mat') % From SCO Estimation DSCC 2016 (Partha, Sergio)
Q = 1000; %mAh  *Acceptable range: 2000 ± 1000 mAh*
th(1,1) = 1/((Q/1000)*3600); %As 
th(1,2) = 0.0042;
th(1,3) = 0.0036; % Thermal time-constant =hA/mCp = 13*5e-4/30e-3*60 = 0.0036
th(1,4) = 0.0427; % 1/mCp  m=30g, Cp=13J/sm2k

% the "th" Parameters are
%th1 = 1/Q,
% th2=1/RC = 1/ (0.1ohm * 2400) = 0.0042
% **Q_battery (A.s) = CV = C*3v = 7200 A.s =>>> C = 2400F  >>> th2 = 1/(0.1*2400) = 0.0042
%th3 = Ah/mCp = 0.0036
%th5 = 1/C;
%th6 = ohmic resistance =~ 0.27 ohm

dt = 5; % [Sec] 
Device_LifeTime = 6.5; % ± 3 [years]
time = 1:dt:Device_LifeTime*3600*24*365; %[sec] ~ 6.5 years of battery operation
%LifeTime_Years = (dt*length(time))/(3600*24*365)% [Years]


% Begin by initiating the state and input values
initSOC = 0.95 %0.95; % 100% SOC (x1) - Fully charged battery
x2 = 0; % when the battery is not initially under the load, there is no voltage drop accross the RC pair in the ECM model
x3 = 37+273.15;
u = -25e-6; % Constant load current ± 25e-6 Amp
u_shock = -3; % Defibrillation current ± 2 Amp

% Shock Pulses Configuration
config = struct();
config(1).shocks = 2; % Number of shocks in the first year
config(1).first_shock_time_sec = 2*30*24*3600; % Time in seconds for the first shock in the first year
config(2).shocks = 2; % Number of shocks in the 2nd year
config(2).first_shock_time_sec = 2*30*24*3600; % Time in seconds for the first shock in the 2nd year
config(3).shocks = 2; % Number of shocks in the 3rd year
config(3).first_shock_time_sec = 2*30*24*3600; % Time in seconds for the first shock in the 3rd year
config(4).shocks = 3; % Number of shocks in the 4th year
config(4).first_shock_time_sec = 2*30*24*3600; % Time in seconds for the first shock in the 4th year
config(5).shocks = 3; % Number of shocks in the 5th year
config(5).first_shock_time_sec = 2*30*24*3600; % Time in seconds for the first shock in the 5th year
config(6).shocks = 3; % Number of shocks in the 6th year
config(6).first_shock_time_sec = 2*30*24*3600; % Time in seconds for the first shock in the 6th year

% Pulse width and interval settings (global or per year)
pulse_width_sec = 10; % Pulse width in seconds
interval_between_shocks_sec = 3*30*24*3600; % Interval between shocks in seconds

% Simulation Parameters
total_time_seconds = Device_LifeTime * 3600 * 24 * 365; % Total simulation time in seconds

% Initialize the current array with value 'u'
current = u * ones(1, total_time_seconds / dt);

% Loop through each year
for year = 1:length(config)
    % Skip if no shocks are defined for the year
    if ~isfield(config(year), 'shocks') || config(year).shocks == 0
        continue;
    end

    % Calculate the start time for the first shock in the year
    year_start_sec = (year - 1) * 3600 * 24 * 365;
    shock_start_sec = year_start_sec + config(year).first_shock_time_sec;

    % Apply the shocks for the year
    for shock_num = 1:config(year).shocks
        shock_start_index = round(shock_start_sec / dt) + 1;

        % Apply the shock for the pulse width duration
        for pulse_sec = 0:dt:pulse_width_sec
            current_index = shock_start_index + round(pulse_sec / dt);
            if pulse_sec < pulse_width_sec
                current(current_index) = u_shock;
            end
        end

        % Update the start time for the next shock
        shock_start_sec = shock_start_sec + pulse_width_sec + interval_between_shocks_sec;
    end
end



% Initialize the three states of the thermal/electrochemical model
x = zeros(3,length(time));
x(:,1) = [initSOC; x2; x3];  % [State of Charge, Voltage Transient, Temp]

% Initialize the output of the thermal/electrochemical model
y = zeros(1,length(time));
y(1) = ocv(100001-(round(x(1,1)*100000)+1));  % Initial Terminal Voltage
% Run a forward Euler solver of the thermal/electrochemical model
for i=1:(length(time)-1)
    x(1,i+1) = x(1,i) + (th(1)*current(i))*dt;
    x(2,i+1) = x(2,i) + (-th(2)*x(2,i)+current(i))*dt;
    x(3,i+1) = x(3,i) + (th(3)*(ambientT-x(3,i))+th(4)*th(6)*current(i)^2+...
               th(4)*(th(7)+th(8)*x(1,i)+th(9)*x(1,i)^2+th(10)*x(1,i)^3+...
               th(11)*x(1,i)^4)*x(3,i)*current(i))*dt;
           

    % Check if SOC (x1) falls below zero
    if x(1,i) < 0.01
        % Truncate arrays to the current time step as simulation stops
        x = x(:, 1:i);
        y = y(1:i);
        current = current(1, 1:i);  % Truncate the current array
        time = time(1:i);
        break; % Exit the loop
    end

    if x(1,i+1) > 0.001 && x(1,i+1) < 0.999 % 0.001<SOC<0.999
        y(i+1) = ocv(100001-(round(x(1,i)*100000)+1))+th(5)*x(2,i+1) + th(6)*current(i+1);
    else
        y(i+1) = ocv(1);
    end
    if y(i+1) >= 3.4
        y(i+1)= 0;
    end

end

T = x(3,:);


subplot(2,1,1)
plot (time/(3600*24*365),y,'LineWidth',1.5,'Color','k');
ylabel('Terminal voltage [v]','FontSize',13,'FontWeight','bold')
xlabel('time [year]','FontSize',12,'FontWeight','bold')
grid on

subplot(2,1,2)
plot (time/(3600*24*365),100-(x(1,:)*100),'LineWidth',1.5,'Color','k');
ylabel('DOD [%]','FontSize',13,'FontWeight','bold')
xlabel('time [year]','FontSize',12,'FontWeight','bold')
grid on


figure (2)
plot (time/(3600*24*365),x(3,:)-273.15,'LineWidth',1.5,'Color','k');
ylabel('Battery temperature [C]','FontSize',13,'FontWeight','bold')
xlabel('time [year]','FontSize',12,'FontWeight','bold')
grid on

figure (3)
plot (time/(3600*24*365),current,'LineWidth',1.5)
ylabel('Battery discharge current [A]','FontSize',12,'FontWeight','bold')
xlabel('time [year]','FontSize',12,'FontWeight','bold')
grid on

figure (4)
plot (soc,ocv,'LineWidth',1.5)
title('Li/SVO-CFx Characteristic curve')
ylabel('OCV [V]','FontSize',12,'FontWeight','bold')
xlabel('SOC [%]','FontSize',12,'FontWeight','bold') 
grid on


%% Saving the results in the Matlab file directory
% Get the directory of the current .m file
currentFolder = fileparts(which(mfilename));

% Name of the folder to save the plots and workspace
saveFolder = 'Saved_Data_Results';

% Full path of the directory to save data
fullSavePath = fullfile(currentFolder, saveFolder);

% Create the directory if it does not exist
if ~exist(fullSavePath, 'dir')
    mkdir(fullSavePath);
end

% Check if the directory was created successfully
if ~exist(fullSavePath, 'dir')
    error('Failed to create the directory: %s', fullSavePath);
end

% Save all open figures
figHandles = findobj('Type', 'figure');
for i = 1:length(figHandles)
    % Check if the figure handle is valid and the figure is open
    if ishandle(figHandles(i))
        figureName = sprintf('Figure%d', i);
        % Save as .png
        pngPath = fullfile(fullSavePath, [figureName '.png']);
        saveas(figHandles(i), pngPath);
    else
        warning('Skipped saving an invalid or closed figure handle.');
    end
end

