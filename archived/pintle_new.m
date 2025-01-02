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
pInj = 60 * 0.95;  
pComb = 35;        
dComb = 76;         
OF = 2.9;          

% Fluid Properties
oxidizer = 'N2O';   
fuel = 'Ethanol';  
dmOx = 1.3502;        
cdOx = 0.65;        
tankTemperature = 298; 

% Fuel Specific Parameters
dmFuel = dmOx/OF;  
cdFuelHole = 0.611;  

% Geometric Constraints
maxChamberToPintleRatio = 5; 
minChamberToPintleRatio = 3;
dPintle = linspace(dComb/maxChamberToPintleRatio, dComb/minChamberToPintleRatio, 5);
rPintle = dPintle/2;

%% Fluid characteristics and areas computation

% Oxidizer density computation
rhoOx = PropSI('D', 'P', (pInj-6)*1e5, 'T', 273, oxidizer); % cold refueling 0°C
aOx = dmOx/(cdOx*sqrt(2*rhoOx*(pInj-pComb)*1e5)); % [m]
aOx = aOx*1e6; % [mm]

% Fuel density computation
rhoFuel = PropSI('D', 'P', (pInj-6)*1e5, 'T', 490, fuel); %assuming 6 bar pressure drop
aFuel = dmFuel/(cdFuelHole*sqrt(2*rhoFuel*(pInj-pComb)*1e5)); % [m]
aFuel = aFuel*1e6; % [mm]

%% Pintle configuration: N2O

rOuterOx = zeros(size(rPintle));

for i = 1:length(rPintle)
    rOuterOx(i) = sqrt( (aOx)/pi + (rPintle(i)*1e-3)^2 );
end

%% Pintle configuration: Ethanol

% Hole diameter options (meters)
dFuelHoleOptions = [0.3, 0.4, 0.5, 0.6, 0.8] * 1e-3;
rowOptions = [1, 2, 3, 4, 5];

% Initialization of variables
total = struct();
vFuel = zeros(length(dFuelHoleOptions), length(rowOptions), length(rPintle));
TMR = zeros(length(dFuelHoleOptions), length(rowOptions), length(rPintle));
alpha = zeros(length(dFuelHoleOptions), length(rowOptions), length(rPintle));
BF = zeros(length(rPintle), length(dFuelHoleOptions));
rOuterOx = zeros(size(rPintle));




% Prepare storage for results
resultsMatrix = cell(length(dFuelHoleOptions), length(rowOptions), length(rPintle));

%rPintle = rPintle*1e-3;
for i = 1:length(dFuelHoleOptions)

    % Compute hole area
    total.holeArea(i) = pi * (dFuelHoleOptions(i)/2)^2;


    % Compute exact number of holes required and then round to the closest integer
    total.exactHoles(i) = dmFuel / (cdFuelHole * total.holeArea(i) * sqrt(2*rhoFuel*(pInj - pComb)*1e5));
    total.maxHoles(i) = round(total.exactHoles(i));

    % Analyze distribution across different row configurations
    for j = 1:length(rowOptions)

        % Compute holes per row (rounded up)
        total.holesPerRow(i, j) = round(total.maxHoles(i) / rowOptions(j));

        % Calculate total holes
        total.totalHoles(i, j) = rowOptions(j) * total.holesPerRow(i,j);

        % Calculate actual mass flow rate
        total.actualMassFlow(i, j) = total.totalHoles(i,j) * total.holeArea(i) * cdFuelHole * sqrt(2*rhoFuel*(pInj - pComb)*1e5);

        % Calculate angular spacing between holes
        total.angularSpacing(i, j) = 360 / total.holesPerRow(i,j);  % Degrees

        % Convert angular spacing to radians (to calculate arcDistance)
        total.angularSpacingRad(i, j) = deg2rad(total.angularSpacing(i,j));

        % Iterate through rPintle values for each row
        for k = 1:length(rPintle)
            % Calculate arc distance for the current rPintle value
            total.arcDistance(i, j, k) = (rPintle(k)*1e-3) * total.angularSpacingRad(i, j)  - dFuelHoleOptions(i);
            
            % External radius of the annulus crown
            total.rOuterOx(k) = sqrt( (aOx)/pi + (rPintle(k)*1e-3)^2 );

            total.voOx(k) = dmOx / (rhoOx * ((total.rOuterOx(k))^2 - (rPintle(k)*1e-3)^2) * pi);
            total.vFuel(i, j) = total.actualMassFlow(i, j)/(rhoFuel*cdFuelHole*total.totalHoles(i, j)*(pi*(dFuelHoleOptions(i)/2)^2)); % * ???

            total.TMR(i, j, k) = (total.actualMassFlow(i, j)*total.vFuel(i, j))/(dmOx*total.voOx(k));
            total.alpha(i, j, k) = rad2deg(acos(1/(1 + total.TMR(i, j, k))));
            total.BF(i, j, k) = (total.totalHoles(i, j)*dFuelHoleOptions(i))/(pi*rPintle(k)*1e-3);
            
            resultsMatrix{i,j,k} = struct(...
                'holeDiameter', dFuelHoleOptions(i) * 1000, ... % mm
                'holeArea', total.holeArea(i) * 1e6, ... % mm²
                'nRows', rowOptions(j), ...
                'holesPerRow', total.holesPerRow(i, j), ...
                'angularSpacing', total.angularSpacing(i, j), ...
                'totalHoles', total.totalHoles(i, j), ...
                'actualMassFlow', total.actualMassFlow(i, j), ...
                'percentError', abs(total.actualMassFlow(i, j) - dmFuel) / dmFuel * 100, ...
                'arcDistance', total.arcDistance(i, j, k), ...
                'rPintle', rPintle(k) ... % Add rPintle value to the structure
            );
        end
    end
end

%% Results

fprintf('-----------------------------------------------------------\n');
fprintf('Fuel Internal with Holes and N2O Annulus \n\n\n');
fprintf('Computed fuel injection area: %.4f mm²\n', aFuel*1e6); 
fprintf('---------------------------------------------------------------------------------------------------------------------------------------------------\n');

fprintf('Fuel Discrete Hole Configuration Analysis\n\n');
fprintf('---------------------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf('| Hole Dia. | Hole Area  | Rows | Holes/Row   | Total Holes  | Mass Flow (kg/s) | Error (%%) || Angular Delta | Arc Distance (mm)  | rPintle (mm)    |\n');
fprintf('|-------------------------------------------------------------------------------------------------------------------------------------------------|\n');
for i = 1:length(dFuelHoleOptions)
    for j = 1:length(rowOptions)
        for k = 1:length(rPintle)
            res = resultsMatrix{i,j,k};

            fprintf('| %.2f mm   | %.4f mm² |   %d  |     %4d    |     %4d     |      %.4f      |   %2.2f%%   ||   %6.2f deg  |     %6.3f mm      |   %6.2f mm   |\n', ...
                res.holeDiameter, ...
                res.holeArea, ...
                res.nRows, ...
                res.holesPerRow, ...
                res.totalHoles, ...
                res.actualMassFlow, ...
                res.percentError, ...
                res.angularSpacing, ...
                res.arcDistance, ...
                res.rPintle ...
            );
        end
    end
    fprintf('|-------------------------------------------------------------------------------------------------------------------------------------------------|\n\n');
end

fprintf('Computed oxidizer injection area: %.4f mm²\n', aOx*1e6);

fprintf('-----------------------------------------------------\n');
fprintf('|   rPintle (mm)   |   rOuterOx (mm)   |\n');
fprintf('|------------------|-------------------|\n');
for i = 1:length(rPintle)
    fprintf('|     %6.2f      |       %6.2f       |\n', rPintle(i), total.rOuterOx(i)*1e3); % Convert to mm
end
fprintf('-----------------------------------------------------\n');

%% Performance

% % Initialize variables to store results

voOx = zeros(length(dFuelHoleOptions), length(rowOptions), length(rPintle)); 
vFuel = zeros(length(dFuelHoleOptions), length(rowOptions), length(rPintle));
TMR = zeros(length(dFuelHoleOptions), length(rowOptions), length(rPintle));
alpha = zeros(length(dFuelHoleOptions), length(rowOptions), length(rPintle));
BF = zeros(length(rPintle), length(dFuelHoleOptions));

for i = 1:length(dFuelHoleOptions)
    voOx(i) = dmOx(i)/(rhoOx*( (doOx(i)*1e-3/2)^2 - (data.ox.diOx(i)*1e-3/2)^2 )*pi);
    vFuel(i) = dmFuelComputed(i)/(rhoFuel*cdFuelHole*nFuel(i)*(pi*(dFuelHole(i)/2)^2));
    TMR(i) = (dmFuelComputed(i)*vFuel(i))/(dmOx*voOx);
    alpha(i) = rad2deg(acos(1/(1+TMR(i))));
    for j=1:length(dInternalOx)
        BF(j,i) = (nFuel(i)*dFuelHole(i))/(pi*dInternalOx(j)*1e-3);
    end
end

for i = 1:length(dFuelHoleOptions)
    for j = 1:length(rowOptions)
        for k = 1:length(rPintle)
            
            % Access actualMassFlow from the structure
            actualMassFlow = total.actualMassFlow(i,j);

            % Compute oxidizer velocity voOx
            voOx(i, j) = dmOx / (rhoOx * ((rOuterOx(k)*1e-3)^2 - (rPintle(k)*1e-3)^2) * pi);
            
            % Compute fuel velocity vFuel
            vFuel(i, j, k) = total.mdot.actualMassFlow(i,j) / ...
                             (rhoFuel * cdFuelHole * result.totalHoles * (pi * (result.holeDiameter*1e-3 / 2)^2));
            
            % Compute Thrust Momentum Ratio (TMR)
            TMR(i, j, k) = (total.mdot.actualMassFlow(i,j) * vFuel(i, j, k)) / (dmOx * voOx(i, j, k));
            
            % Compute spray angle alpha
            alpha(i, j, k) = rad2deg(acos(1 / (1 + TMR(i, j, k))));

            % Print results for current configuration
            fprintf('Configuration (dFuelHoleOption=%.3f mm, rowOption=%d, rPintle=%.3f mm):\n', ...
                dFuelHoleOptions(i)*1e3, rowOptions(j), rPintle(k)*1e3);
            fprintf('  voOx = %.4f m/s\n', voOx(i, j, k));
            fprintf('  vFuel = %.4f m/s\n', vFuel(i, j, k));
            fprintf('  TMR = %.4f\n', TMR(i, j, k));
            fprintf('  alpha = %.2f degrees\n\n', alpha(i, j, k));
        end
    end

    % Compute blockage factor (BF) for each `dInternalOx`
    for j = 1:length(rPintle)
        % Compute nFuel using `result.totalHoles`
        nFuel = result.totalHoles; % This uses the number of holes from resultsMatrix
        
        % Compute blockage factor BF
        BF(j, i) = (nFuel * dFuelHoleOptions(i)) / (pi * rPintle(j) * 1e-3);

        % Print blockage factor
        fprintf('Blockage Factor (dFuelHoleOption=%.3f mm, dInternalOx=%.3f mm):\n', ...
            dFuelHoleOptions(i)*1e3, rPintle(j)*1e3);
        fprintf('  BF = %.4f\n\n', BF(j, i));
    end
end


% fprintf('\n------------------------DATA RECORD for DESIGN-------------------------\n\n');
fprintf('OXIDIZER \n');
fprintf('The oxidizer area of the annulus is %.2f mm^2 with the following diameter: \n', aOx)
for i = 1:length(dInternalOx)
    fprintf(' Inner diameter of %.2f mm and an outer diameter of %.2f mm \n', dInternalOx(i), doOx(i))
end
fprintf('The minimum pipe diameter for the oxidizer is %.2f mm \n\n', minPipeDiameterOx);

fprintf('FUEL \n');
fprintf('The number of holes for a desired area of %.2f mm^2 the fuel is: \n ', aFuel)
for i = 1:length(dFuelHole)
    fprintf(' %d with a diameter of %.2f mm for a mass flow rate of %.4f kg/s \n ', nFuel(i), dFuelHole(i)*1e3, dmFuelComputed(i))
end
fprintf('The minimum pipe diameter for the fuel is %.2f mm \n\n', minPipeDiameterFuel);

fprintf('\n------------------------PERFORMANCE-------------------------\n\n');
fprintf('For what concern the performance of the pintle injector, what was found out is that TMR depends on fuel mass flow rate: \n')
for i = 1:length(dFuelHole)
    fprintf(' TMR = %.4f with a mass flow rate of %.4f kg/s \n', TMR(i), dmFuelComputed(i))
end
fprintf('\nWhile the respective spray angles are: \n')
for i = 1:length(dFuelHole)
    fprintf(' alpha = %.2f° for TMR %.4f \n', alpha(i), TMR(i))
end

%% Function

function result = PropSI(output, input1, value1, input2, value2, fluid)
    result = py.CoolProp.CoolProp.PropsSI(output, input1, value1, input2, value2, fluid);
end
