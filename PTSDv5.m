clc; clear

% DESCRIPTION:
%   Script for analyzing and designing pintle injector configurations 
%   with two different geometrical arrangements (Oxidizer Internal/Fuel External 
%   and Fuel Central/Oxidizer External)
%
% --------------------------------------------------------------------------------------
% NOMENCLATURE
% pInj             : Injection pressure, considering a 5% pressure drop in lines [bar]
% pComb            : Combustion chamber pressure [bar]
% dComb            : Combustion chamber diameter [mm]
% OF               : Oxidizer-to-Fuel ratio (dimensionless) [-]
% dmOx             : Mass flow rate of oxidizer [kg/s]
% cdOx             : Discharge coefficient for oxidizer injector [-]
% tankTemperature  : Oxidizer tank temperature [K]
% dmFuel           : Mass flow rate of fuel [kg/s]
% cdFuelHole       : Discharge coefficient for fuel injector holes [-]

% --------------------------------------------------------------------------------------
%
% OUTPUT:
%   Comparative analysis of pintle injector configurations with various 
%   geometric and performance parameters
%
% EXAMPLE:
%   Run the script directly to generate comparative output
%
% CALLED FUNCTIONS:
%   - PropSI (custom function interfacing with CoolProp)
%   - fChockParamFun.mat
%
% ASSUMPTIONS AND LIMITATIONS:
%   - Uses CoolProp for thermodynamic property calculations
%   - Assumes specific injection and combustion conditions
%   - Uses simplified NHNE model for injection characterization
%
% REFERENCES:
%   - CoolProp thermodynamic library
%   - NHNE injection model
%   - SPI injection model
%   - Naming scheme as reported on Notion dedicated page

% REVISION HISTORY:
%   29/11/2024 - Christian Ferracane
%       * Added script header and reorganized constants
%   26/11/2024 - Christian Ferracane
%       * Added the two configurations, implemented PropSI (scrapped broken
%       interpolation).
%   25/11/2024 - Chiara Conte
%       * First generation of the code from previous year injector script
%       and plots.
%
% LICENSE:
%   Copyright © 2023, Skyward Experimental Rocketry, PRP department
%   All rights reserved
%   SPDX-License-Identifier: GPL-3.0-or-later
%
% CONTACT:
%   Skyward Experimental Rocketry: info@skywarder.eu

%% Parameters definition:

% Injection and Combustion Parameters
pInj_fuel = 60 * 0.95 - 6;  % Already with the pressure drop
pInj_ox = 50 * 0.95;
pComb = 35;        
dComb = 76;         
OF = 2.9;          

% Fluid Properties
oxidizer = 'N2O';   
fuel = 'Ethanol';  
dmOx = 1.3502;        
cdAnnulus = 0.65;        
tankTemperature = 298; 
oxTemp = 273.15 + 0;

% Fuel Specific Parameters
dmFuel = dmOx/OF;  
cdHole = 0.611;  

% Geometric Constraints
maxChamberToPintleRatio = 5; 
minChamberToPintleRatio = 3;
dPintle = linspace(dComb/maxChamberToPintleRatio, dComb/minChamberToPintleRatio, 5);
rPintle = dPintle/2;

%% Fluid characteristics and areas computation

% Oxidizer density computation
rhoOx = PropSI('D', 'P', pInj_ox*1e5, 'T', oxTemp, oxidizer); % TODO: determine if cold refueling is problematic
aOx = dmOx/(cdAnnulus*sqrt(2*rhoOx*(pInj_ox-pComb)*1e5)); % [m]
% aOx = dmOx/(cdAnnulus*sqrt(2*rhoOx*(20)*1e5)); % [m]
%aOx = aOx*1e6; % [mm]

% Fuel density computation
rhoFuel = PropSI('D', 'P', (pInj_fuel-6)*1e5, 'T', 490, fuel); %assuming 6 bar pressure drop
aFuel = dmFuel/(cdHole*sqrt(2*rhoFuel*(pInj_fuel-pComb)*1e5)); % [m]
% aFuel = aFuel*1e6; % [mm]

%% Pintle configuration and performance Fuel centrale e Ox esterno

% Hole diameter options (meters)
dFuelHoleOptions = [0.5, 0.6, 0.8, 1, 1.2] * 1e-3;
rowOptions = [1, 2, 3, 4, 5];

% Variables initialization
holeArea = zeros(size(dFuelHoleOptions));
nHoles = zeros(size(dFuelHoleOptions));
maxnHoles = zeros(size(dFuelHoleOptions));

rOuterOx = zeros(size(rPintle));
voOx = zeros(size(rPintle));


holesPerRows = zeros(length(dFuelHoleOptions), length(rowOptions));
totalHoles = zeros(length(dFuelHoleOptions), length(rowOptions));
actualDm = zeros(length(dFuelHoleOptions), length(rowOptions));
angularSpacing = zeros(length(dFuelHoleOptions), length(rowOptions));
angularSpacingRad = zeros(length(dFuelHoleOptions), length(rowOptions));
vFuel = zeros(length(dFuelHoleOptions), length(rowOptions));
BF = zeros(length(dFuelHoleOptions), length(rowOptions));
percError = zeros(length(dFuelHoleOptions), length(rowOptions));

arcDistance = zeros(length(dFuelHoleOptions), length(rowOptions), length(rPintle));
TMR = zeros(length(dFuelHoleOptions), length(rowOptions), length(rPintle));
alpha = zeros(length(dFuelHoleOptions), length(rowOptions), length(rPintle));

% Prepare storage for results
resultsMatrix = cell(length(dFuelHoleOptions), length(rowOptions), length(rPintle));

for i = 1:length(dFuelHoleOptions)

    holeArea(i) = pi*(dFuelHoleOptions(i)/2)^2;
    nHoles(i) = dmFuel / (cdHole * holeArea(i) * sqrt(2*rhoFuel*(pInj_fuel - pComb)*1e5) );
    maxnHoles(i) = floor(nHoles(i));

    for j = 1:length(rowOptions)

        holesPerRows(i, j) = floor(nHoles(i) / rowOptions(j));
        totalHoles(i, j) = holesPerRows(i, j) * rowOptions(j);
        actualDm(i, j) = cdHole * totalHoles(i, j) * holeArea(i) * sqrt(2*rhoFuel*(pInj_fuel - pComb)*1e5);
        percError(i, j) = abs(actualDm(i, j) - dmFuel) / dmFuel * 100;
        angularSpacing(i, j) = 360 / holesPerRows(i, j);
        angularSpacingRad(i, j) = deg2rad(angularSpacing(i,j));

        for k = 1:length(rPintle)
            arcDistance(i, j, k) = (((rPintle(k)*1e-3) * angularSpacingRad(i, j)) - dFuelHoleOptions(i))*1e3; %mm

            rOuterOx(k) = sqrt( (aOx)/pi + (rPintle(k)*1e-3)^2 );
            %voOx(k) = dmOx / ( rhoOx * (( (rOuterOx(k))^2 - (rPintle(k)*10^(-3))^2 ) * pi ));
            voOx(k) = dmOx / (rhoOx * pi* ((rOuterOx(k))^2 - (rPintle(k)*1e-3)^2));
            vFuel(i, j) = actualDm(i, j)/(rhoFuel * cdHole * totalHoles(i, j) * holeArea(i));
            TMR(i, j, k) = (actualDm(i, j) * vFuel(i, j) )/(dmOx*voOx(k));

            alpha(i, j, k) = rad2deg(acos(1/(1 + TMR(i, j, k))));
            BF(i, j, k) = (holesPerRows(i, j) * dFuelHoleOptions(i)) / (2 * pi * rPintle(k)*1e-3);

            resultsMatrix{i,j,k} = struct(...
                'holeDiameter', dFuelHoleOptions(i) * 1000, ... % mm
                'holeArea', holeArea(i) * 1e6, ... % mm²
                'rowOptions', rowOptions(j), ...
                'holesPerRows',  holesPerRows(i, j), ...
                'angularSpacing', angularSpacing(i, j), ...
                'totalHoles', totalHoles(i, j), ...
                'actualMassFlow', actualDm(i, j), ...
                'percentError', percError(i, j), ...
                'arcDistance', arcDistance(i, j, k), ...
                'rPintle', rPintle(k), ... 
                'annulusWidth', (rOuterOx(k) * 1e3) - rPintle(k), ...
                'voOx', voOx(k), ...
                'vFuel', vFuel(i,j), ...
                'TMR', TMR(i, j, k), ...
                'alpha', alpha(i, j, k), ...
                'BF', BF(i, j, k) ...
            );
        end
    end
end

%% Plots design parameters

fprintf('-----------------------DESIGN PARAMETERS-------------------------------\n');
fprintf('-----------------------------------------------------------\n');
fprintf('Fuel Internal with Holes and N2O Annulus \n\n\n');
fprintf('Computed fuel injection area: %.4f mm²\n', aFuel*1e6); 
fprintf('Computer oxidizer injection area: %.4f mm²\n',aOx*1e6)
fprintf('-----------------------------------------------------------\n');

fprintf('Fuel Discrete Hole Configuration Analysis\n\n');
fprintf('------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf('| Hole Dia. | Hole Area  | Rows | Holes/Row   | Total Holes  | Mass Flow (kg/s) | Error (%%) || Angular Delta | Arc Distance (mm)  | rPintle (mm) | Annulus Width (mm)  |\n');
fprintf('|-----------------------------------------------------------------------------------------------------------------------------------------------------------------|\n');

for i = 1:length(dFuelHoleOptions)
    for j = 1:length(rowOptions)
        for k = 1:length(rPintle)

            fprintf('| %.2f mm   | %.4f mm² |   %d  |     %4d    |     %4d     |      %.4f      |   %2.2f%%   ||   %6.2f deg  |     %6.3f mm      |   %6.2f mm   | %6.2f mm \n', ...
               resultsMatrix{i,j,k}.holeDiameter, ...
               resultsMatrix{i,j,k}.holeArea, ...
               resultsMatrix{i,j,k}.rowOptions, ...
               resultsMatrix{i,j,k}.holesPerRows, ...
               resultsMatrix{i,j,k}.totalHoles, ...
               resultsMatrix{i,j,k}.actualMassFlow, ...
               resultsMatrix{i,j,k}.percentError, ...
               resultsMatrix{i,j,k}.angularSpacing, ...
               resultsMatrix{i,j,k}.arcDistance, ...
               resultsMatrix{i,j,k}.rPintle, ...
               resultsMatrix{i,j,k}.annulusWidth ...
            );
        end
    end
    fprintf('|------------------------------------------------------------------------------------------------------------------------------------------------------------------|\n\n\n\n\n');
end


%% Plot Performances

fprintf('--------------------------------------------------PERFORMANCE PARAMETERS--------------------------------------------------\n');
fprintf('| Hole Dia. | Hole Area  | rPintle (mm)  |   v_oxid (m/s)  |   v_fuel (m/s) |    TMR     |   alpha   |     BF     | Rows |  Holes/Row  | Total Holes |  \n');
fprintf('|--------------------------------------------------------------------------------------------------------------------------------------------------------|\n');
for i = 1:length(dFuelHoleOptions)
    for j = 1:length(rowOptions)
        for k = 1:length(rPintle)

            fprintf('| %.2f mm   | %.4f mm² |   %6.2f mm   |   %7.4f m/s   |   %7.4f m/s  |   %6.4f   |   %5.2f°  |   %6.4f   |   %d  |    %4d     |    %4d    |   \n', ...
               resultsMatrix{i,j,k}.holeDiameter, ...
               resultsMatrix{i,j,k}.holeArea, ...
               resultsMatrix{i,j,k}.rPintle, ...
               resultsMatrix{i,j,k}.voOx, ...
               resultsMatrix{i,j,k}.vFuel, ...
               resultsMatrix{i,j,k}.TMR, ...
               resultsMatrix{i,j,k}.alpha, ...
               resultsMatrix{i,j,k}.BF, ...
               resultsMatrix{i,j,k}.rowOptions, ...
               resultsMatrix{i,j,k}.holesPerRows, ...
               resultsMatrix{i,j,k}.totalHoles ...
            );
        end
    end
    fprintf('|-------------------------------------------------------------------------------------------------------------------------------------------------|\n\n');
end

%% Excel file creation
% File name for Excel
filename = 'PintleInjectorDesign.xlsx';

% PropSI Results
rhoOx = PropSI('D', 'P', pInj_ox*1e5, 'T', oxTemp, oxidizer); 
rhoFuel = PropSI('D', 'P', (pInj_fuel-6)*1e5, 'T', 490, fuel); 

% Write Design Parameters
designParams = {'Computed Fuel Injection Area (mm²)', aFuel*1e6; ...
                'Computed Oxidizer Injection Area (mm²)', aOx*1e6; ...
                'Oxidizer Density (kg/m³)', rhoOx; ...
                'Fuel Density (kg/m³)', rhoFuel};
writecell({'DESIGN PARAMETERS'}, filename, 'Sheet', 'Design', 'Range', 'A1');
writecell(designParams, filename, 'Sheet', 'Design', 'Range', 'A3');

fuelDiscreteAnalysis = {'Hole Dia. (mm)', 'Hole Area (mm²)', 'Rows', ...
                        'Holes/Row', 'Total Holes', 'Mass Flow (kg/s)', ...
                        'Error (%)', 'Angular Delta (deg)', ...
                        'Arc Distance (mm)', 'rPintle (mm)', 'Annulus Width (mm)'};

for i = 1:length(dFuelHoleOptions)
    for j = 1:length(rowOptions)
        for k = 1:length(rPintle)
            % Convert numeric array into a cell array before concatenating
            dataRow = {resultsMatrix{i,j,k}.holeDiameter, ...
                       resultsMatrix{i,j,k}.holeArea, ...
                       resultsMatrix{i,j,k}.rowOptions, ...
                       resultsMatrix{i,j,k}.holesPerRows, ...
                       resultsMatrix{i,j,k}.totalHoles, ...
                       resultsMatrix{i,j,k}.actualMassFlow, ...
                       resultsMatrix{i,j,k}.percentError, ...
                       resultsMatrix{i,j,k}.angularSpacing, ...
                       resultsMatrix{i,j,k}.arcDistance, ...
                       resultsMatrix{i,j,k}.rPintle, ...
                       resultsMatrix{i,j,k}.annulusWidth};
            fuelDiscreteAnalysis = [fuelDiscreteAnalysis; dataRow];
        end
    end
end

writecell({'FUEL DISCRETE HOLE CONFIGURATION ANALYSIS'}, filename, 'Sheet', 'Fuel Analysis', 'Range', 'A1');
writecell(fuelDiscreteAnalysis, filename, 'Sheet', 'Fuel Analysis', 'Range', 'A3')

% Write Performance Parameters
performanceParams = {'Hole Dia. (mm)', 'Hole Area (mm²)', 'rPintle (mm)', ...
                     'v_oxid (m/s)', 'v_fuel (m/s)', 'TMR', 'alpha (deg)', ...
                     'BF', 'Rows', 'Holes/Row', 'Total Holes'};

for i = 1:length(dFuelHoleOptions)
    for j = 1:length(rowOptions)
        for k = 1:length(rPintle)
            % Convert numeric array into a cell array
            dataRow = {resultsMatrix{i,j,k}.holeDiameter, ...
                       resultsMatrix{i,j,k}.holeArea, ...
                       resultsMatrix{i,j,k}.rPintle, ...
                       resultsMatrix{i,j,k}.voOx, ...
                       resultsMatrix{i,j,k}.vFuel, ...
                       resultsMatrix{i,j,k}.TMR, ...
                       resultsMatrix{i,j,k}.alpha, ...
                       resultsMatrix{i,j,k}.BF, ...
                       resultsMatrix{i,j,k}.rowOptions, ...
                       resultsMatrix{i,j,k}.holesPerRows, ...
                       resultsMatrix{i,j,k}.totalHoles};
            performanceParams = [performanceParams; dataRow];
        end
    end
end

writecell({'PERFORMANCE PARAMETERS'}, filename, 'Sheet', 'Performance', 'Range', 'A1');
writecell(performanceParams, filename, 'Sheet', 'Performance', 'Range', 'A3');

%% Pintle configuration and performance Fuel esterno e Ox centrale

% Hole diameter options (meters)
dOxHoleOptions = [0.5, 0.6, 0.8, 1, 1.2] * 1e-3;
rowOptions = [1, 2, 3, 4, 5];

% Variables initialization
holeArea = zeros(size(dOxHoleOptions));
nHoles = zeros(size(dOxHoleOptions));
maxnHoles = zeros(size(dOxHoleOptions));

rOuterFuel = zeros(size(rPintle));
vFuel = zeros(size(rPintle));


holesPerRows = zeros(length(dOxHoleOptions), length(rowOptions));
totalHoles = zeros(length(dOxHoleOptions), length(rowOptions));
actualDm = zeros(length(dOxHoleOptions), length(rowOptions));
angularSpacing = zeros(length(dOxHoleOptions), length(rowOptions));
angularSpacingRad = zeros(length(dOxHoleOptions), length(rowOptions));
voOx = zeros(length(dOxHoleOptions), length(rowOptions));
BF = zeros(length(dOxHoleOptions), length(rowOptions));
percError = zeros(length(dOxHoleOptions), length(rowOptions));

arcDistance = zeros(length(dOxHoleOptions), length(rowOptions), length(rPintle));
TMR = zeros(length(dOxHoleOptions), length(rowOptions), length(rPintle));
alpha = zeros(length(dOxHoleOptions), length(rowOptions), length(rPintle));

% Prepare storage for results
resultsMatrix = cell(length(dOxHoleOptions), length(rowOptions), length(rPintle));

for i = 1:length(dOxHoleOptions)

    holeArea(i) = pi*(dOxHoleOptions(i)/2)^2;
    nHoles(i) = dmOx / (cdHole * holeArea(i) * sqrt(2*rhoOx*(pInj_ox - pComb)*1e5) );
    maxnHoles(i) = floor(nHoles(i));

    for j = 1:length(rowOptions)

        holesPerRows(i, j) = floor(nHoles(i) / rowOptions(j));
        totalHoles(i, j) = holesPerRows(i, j) * rowOptions(j);
        actualDm(i, j) = cdHole * totalHoles(i, j) * holeArea(i) * sqrt(2*rhoOx*(pInj_ox - pComb)*1e5);
        percError(i, j) = abs(actualDm(i, j) - dmOx) / dmOx * 100;
        angularSpacing(i, j) = 360 / holesPerRows(i, j);
        angularSpacingRad(i, j) = deg2rad(angularSpacing(i,j));

        for k = 1:length(rPintle)
            arcDistance(i, j, k) = (((rPintle(k)*1e-3) * angularSpacingRad(i, j))- dOxHoleOptions(i))*1e3; %mm

            rOuterFuel(k) = sqrt( (aFuel)/pi + (rPintle(k)*1e-3)^2 );
            %voOx(k) = dmOx / ( rhoOx * (( (rOuterOx(k))^2 - (rPintle(k)*10^(-3))^2 ) * pi ));
            vFuel(k) = dmFuel / (rhoFuel * pi* ((rOuterFuel(k))^2 - (rPintle(k)*1e-3)^2));
            voOx(i, j) = actualDm(i, j)/(rhoOx * cdHole * totalHoles(i, j) * holeArea(i));
            TMR(i, j, k) = (actualDm(i, j) * voOx(i, j) )/(dmFuel*vFuel(k));

            alpha(i, j, k) = rad2deg(acos(1/(1 + TMR(i, j, k))));
            BF(i, j, k) = (holesPerRows(i, j) * dOxHoleOptions(i)) / (2 * pi * rPintle(k)*1e-3);

            resultsMatrix{i,j,k} = struct(...
                'holeDiameter', dOxHoleOptions(i) * 1000, ... % mm
                'holeArea', holeArea(i) * 1e6, ... % mm²
                'rowOptions', rowOptions(j), ...
                'holesPerRows',  holesPerRows(i, j), ...
                'angularSpacing', angularSpacing(i, j), ...
                'totalHoles', totalHoles(i, j), ...
                'actualMassFlow', actualDm(i, j), ...
                'percentError', percError(i, j), ...
                'arcDistance', arcDistance(i, j, k), ...
                'rPintle', rPintle(k), ... 
                'annulusWidth', (rOuterFuel(k) * 1e3) - rPintle(k), ...
                'vFuel', vFuel(k), ...
                'voOx', voOx(i,j), ...
                'TMR', TMR(i, j, k), ...
                'alpha', alpha(i, j, k), ...
                'BF', BF(i, j, k) ...
            );
        end
    end
end

%% Plots design parameters

fprintf('-----------------------DESIGN PARAMETERS-------------------------------\n');
fprintf('-----------------------------------------------------------\n');
fprintf('N2O Internal with Holes and fuel Annulus \n\n\n');
fprintf('Computed fuel injection area: %.4f mm²\n', aFuel*1e6); 
fprintf('Computer oxidizer injection area: %.4f mm²\n', aOx*1e6)
fprintf('-----------------------------------------------------------\n');

fprintf('Fuel Discrete Hole Configuration Analysis\n\n');
fprintf('------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf('| Hole Dia. | Hole Area  | Rows | Holes/Row   | Total Holes  | Mass Flow (kg/s) | Error (%%) || Angular Delta | Arc Distance (mm)  | rPintle (mm) | Annulus Width (mm)  |\n');
fprintf('|-----------------------------------------------------------------------------------------------------------------------------------------------------------------|\n');

for i = 1:length(dOxHoleOptions)
    for j = 1:length(rowOptions)
        for k = 1:length(rPintle)

            fprintf('| %.2f mm   | %.4f mm² |   %d  |     %4d    |     %4d     |      %.4f      |   %2.2f%%   ||   %6.2f deg  |     %6.3f mm      |   %6.2f mm   | %6.2f mm \n', ...
               resultsMatrix{i,j,k}.holeDiameter, ...
               resultsMatrix{i,j,k}.holeArea, ...
               resultsMatrix{i,j,k}.rowOptions, ...
               resultsMatrix{i,j,k}.holesPerRows, ...
               resultsMatrix{i,j,k}.totalHoles, ...
               resultsMatrix{i,j,k}.actualMassFlow, ...
               resultsMatrix{i,j,k}.percentError, ...
               resultsMatrix{i,j,k}.angularSpacing, ...
               resultsMatrix{i,j,k}.arcDistance, ...
               resultsMatrix{i,j,k}.rPintle, ...
               resultsMatrix{i,j,k}.annulusWidth ...
            );
        end
    end
    fprintf('|------------------------------------------------------------------------------------------------------------------------------------------------------------------|\n\n\n\n\n');
end



%% Plot Performances

fprintf('--------------------------------------------------PERFORMANCE PARAMETERS--------------------------------------------------\n');
fprintf('| Hole Dia. | Hole Area  | rPintle (mm)  |   v_oxid (m/s)  |   v_fuel (m/s) |    TMR     |   alpha   |     BF     | Rows |  Holes/Row  | Total Holes |  \n');
fprintf('|--------------------------------------------------------------------------------------------------------------------------------------------------------|\n');
for i = 1:length(dFuelHoleOptions)
    for j = 1:length(rowOptions)
        for k = 1:length(rPintle)

            fprintf('| %.2f mm   | %.4f mm² |   %6.2f mm   |   %7.4f m/s   |   %7.4f m/s  |   %6.4f   |   %5.2f°  |   %6.4f   |   %d  |    %4d     |    %4d    |   \n', ...
               resultsMatrix{i,j,k}.holeDiameter, ...
               resultsMatrix{i,j,k}.holeArea, ...
               resultsMatrix{i,j,k}.rPintle, ...
               resultsMatrix{i,j,k}.voOx, ...
               resultsMatrix{i,j,k}.vFuel, ...
               resultsMatrix{i,j,k}.TMR, ...
               resultsMatrix{i,j,k}.alpha, ...
               resultsMatrix{i,j,k}.BF, ...
               resultsMatrix{i,j,k}.rowOptions, ...
               resultsMatrix{i,j,k}.holesPerRows, ...
               resultsMatrix{i,j,k}.totalHoles ...
            );
        end
    end
    fprintf('|-------------------------------------------------------------------------------------------------------------------------------------------------|\n\n');
end

%%

function result = PropSI(output, input1, value1, input2, value2, fluid)
    result = py.CoolProp.CoolProp.PropsSI(output, input1, value1, input2, value2, fluid);
end