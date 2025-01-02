clc; clear;

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

% PARAMETERS DEFINITION

% Injection and Combustion Parameters
pInj = 60 * 0.95;  % Injector pressure after 5% drop [bar]
pComb = 35;        % Combustion chamber pressure [bar]
dComb = 76;        % Combustion chamber diameter [mm]
OF = 2.9;          % Oxidizer-to-fuel ratio (dimensionless)

% Fluid Properties
oxidizer = 'N2O';  % Oxidizer
fuel = 'Ethanol';  % Fuel
dmOx = 1.3502;     % Oxidizer mass flow rate [kg/s]
cdOx = 0.65;       % Discharge coefficient for oxidizer injector
tankTemperature = 298; % Tank temperature [K]

% Fuel Specific Parameters
dmFuel = dmOx / OF;  % Fuel mass flow rate [kg/s]
cdFuelHole = 0.611;  % Discharge coefficient for fuel injector holes

% Geometric Constraints
maxChamberToPintleRatio = 5; 
minChamberToPintleRatio = 3;
dPintle = linspace(dComb / maxChamberToPintleRatio, dComb / minChamberToPintleRatio, 5); % [mm]
rPintle = dPintle / 2; % [mm]

%% FLUID CHARACTERISTICS AND AREA COMPUTATION

% Oxidizer density computation
rhoOx = PropSI('D', 'P', (pInj - 6) * 1e5, 'T', 273, oxidizer); % Cold refueling at 0°C
aOx = dmOx / (cdOx * sqrt(2 * rhoOx * (pInj - pComb) * 1e5)); % [m²]

% Fuel density computation
rhoFuel = PropSI('D', 'P', (pInj - 6) * 1e5, 'T', 490, fuel); % Assuming 6-bar pressure drop
aFuel = dmFuel / (cdFuelHole * sqrt(2 * rhoFuel * (pInj - pComb) * 1e5)); % [m²]

%% PINTLE CONFIGURATION AND PERFORMANCE

% Hole diameter options (meters)
dFuelHoleOptions = [0.3, 0.4, 0.5, 0.6, 0.8] * 1e-3;
rowOptions = [1, 2, 3, 4, 5];

% Initialize storage variables
resultsMatrix = cell(length(dFuelHoleOptions), length(rowOptions), length(rPintle));

for i = 1:length(dFuelHoleOptions)
    holeArea = pi * (dFuelHoleOptions(i) / 2)^2; % Hole area [m²]
    nHoles = dmFuel / (cdFuelHole * holeArea * sqrt(2 * rhoFuel * (pInj - pComb) * 1e5)); 
    maxnHoles = floor(nHoles); % Maximum number of holes

    for j = 1:length(rowOptions)
        holesPerRow = max(1, floor(nHoles / rowOptions(j))); % Holes per row
        totalHoles = holesPerRow * rowOptions(j); % Total number of holes
        actualDmFuel = cdFuelHole * totalHoles * holeArea * sqrt(2 * rhoFuel * (pInj - pComb) * 1e5);
        percentError = abs(actualDmFuel - dmFuel) / dmFuel * 100;
        angularSpacing = 360 / holesPerRow; % Angular delta [degrees]
        angularSpacingRad = deg2rad(angularSpacing); % Angular delta [radians]

        for k = 1:length(rPintle)
            % Arc distance calculation
            arcDistance = rPintle(k) * angularSpacingRad; % [mm]
            
            % Annular dimensions
            rOuterOx = sqrt((aOx * 1e6) / pi + (rPintle(k))^2); % Outer oxidizer radius [mm]
            annulusWidth = rOuterOx - rPintle(k); % Annulus width [mm]
            
            % Fluid velocities
            voOx = dmOx / (rhoOx * pi * ((rOuterOx * 1e-3)^2 - (rPintle(k) * 1e-3)^2)); % Oxidizer velocity [m/s]
            vFuel = actualDmFuel / (rhoFuel * totalHoles * holeArea * 1e6); % Fuel velocity [m/s]
            
            % Thrust Momentum Ratio (TMR)
            TMR = (actualDmFuel * vFuel) / (dmOx * voOx);
            alpha = rad2deg(acos(1 / (1 + TMR))); % Spray angle [degrees]
            
            % Blockage Factor (BF)
            BF = (holesPerRow * dFuelHoleOptions(i) * 1e3) / (2 * pi * rPintle(k)); 

            % Store results
            resultsMatrix{i, j, k} = struct(...
                'holeDiameter', dFuelHoleOptions(i) * 1000, ... % mm
                'holeArea', holeArea * 1e6, ... % mm²
                'rowOptions', rowOptions(j), ...
                'holesPerRow', holesPerRow, ...
                'angularSpacing', angularSpacing, ...
                'arcDistance', arcDistance, ...
                'rPintle', rPintle(k), ...
                'annulusWidth', annulusWidth, ...
                'voOx', voOx, ...
                'vFuel', vFuel, ...
                'TMR', TMR, ...
                'alpha', alpha, ...
                'BF', BF, ...
                'percentError', percentError ...
            );
        end
    end
end

%% FUNCTION DEFINITIONS

function result = PropSI(output, input1, value1, input2, value2, fluid)
    % PropSI wrapper for CoolProp
    result = py.CoolProp.CoolProp.PropsSI(output, input1, value1, input2, value2, fluid);
end
