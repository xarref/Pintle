clc; clear
%% --- Input Parameters ---
% --- Operating Conditions ---
pInj_fuel = 60 * 0.95 - 6;      % Fuel injection pressure [bar]
pInj_ox = 50 * 0.95;            % Oxidizer injection pressure [bar]
pComb = 35;                     % Chamber pressure [bar]
dComb = 76;                     % Chamber diameter [mm]
OF = 2.9;                       % Oxidizer to fuel ratio [-]
% --- Mass Flow Parameters ---
dmOx = 1.3502;                 % Oxidizer mass flow [kg/s]
dmFuel = dmOx / OF;             % Fuel mass flow [kg/s]
% --- Fluid Properties ---
oxidizer = 'N2O';              % Oxidizer fluid (CoolProp name)
fuel = 'Ethanol';              % Fuel fluid (CoolProp name)
tankTemperature = 298;         % Tank temperature [K]
oxTemp = 273.15;               % Oxidizer temperature [K]
% --- Discharge Coefficients ---
cdAnnulus = 0.65;              % Annular injection Cd [-]
cdHole = 0.611;                % Discrete holes Cd [-]
% --- Geometric Parameters ---
maxChamberToPintleRatio = 5;
minChamberToPintleRatio = 3;
dPintle = linspace(dComb / maxChamberToPintleRatio, dComb / minChamberToPintleRatio, 5);
rPintle = dPintle / 2;
dHoleOptions = [0.5, 0.6, 0.8, 1, 1.2] * 1e-3;  % [m] - Hole diameter options
rowOptions = 1:5;                             % [-] - Number of rows options
%% --- Fluid Properties Calculation ---
% --- Get fluid densities using CoolProp ---
try
    rhoOx = PropSI('D', 'P', pInj_ox * 1e5, 'T', oxTemp, oxidizer);
    rhoFuel = PropSI('D', 'P', (pInj_fuel - 6) * 1e5, 'T', 490, fuel); % Assuming fuel temp of 490K for density calculation
catch
    error('CoolProp fluid property calculation failed. Ensure CoolProp is installed and fluids are correctly specified.');
end
% --- Calculate required injection areas ---
aOx = dmOx / (cdAnnulus * sqrt(2 * rhoOx * (pInj_ox - pComb) * 1e5));
aFuel = dmFuel / (cdHole * sqrt(2 * rhoFuel * (pInj_fuel - pComb) * 1e5));
%% --- Initialize Storage ---
configs = {'Fuel_Internal', 'Oxidizer_Internal'};
nConfigs = length(configs);
nHoles = length(dHoleOptions);
nRows = length(rowOptions);
nPintles = length(rPintle);
results = struct(); % Structure to store all results
for c = 1:nConfigs
    results.(configs{c}) = cell(nHoles, nRows, nPintles);
end
%% --- Configuration Calculation Loop ---
for configType = 1:nConfigs
    isFuelInternal = (configType == 1);
    % --- Set fluid properties based on configuration ---
    if isFuelInternal
        fluid = 'Fuel';
        rho = rhoFuel;
        dm = dmFuel;
        pInj = pInj_fuel;
        otherFluidDm = dmOx;
        otherFluidRho = rhoOx;
    else
        fluid = 'Oxidizer';
        rho = rhoOx;
        dm = dmOx;
        pInj = pInj_ox;
        otherFluidDm = dmFuel;
        otherFluidRho = rhoOx;
    end
    % --- Loop through geometric parameters ---
    for i = 1:nHoles
        dHole = dHoleOptions(i);
        holeArea = pi * (dHole / 2)^2;
        for j = 1:nRows
            numRows = rowOptions(j); % Renamed to avoid shadowing nRows (loop index)
            holesPerRow = floor(dm / (cdHole * holeArea * sqrt(2 * rho * (pInj - pComb) * 1e5) * numRows));
            totalHoles = holesPerRow * numRows;
            actualDm = cdHole * totalHoles * holeArea * sqrt(2 * rho * (pInj - pComb) * 1e5);
            angularSpacing = 360 / holesPerRow;
            for k = 1:nPintles
                % --- Calculate geometrical parameters ---
                rPintle_mm = rPintle(k); % Pintle radius in mm for readability
                rOut = sqrt(isFuelInternal * aFuel / pi + (~isFuelInternal * aOx / pi) + (rPintle_mm * 1e-3)^2);
                arcDistance = ((rPintle_mm * 1e-3 * deg2rad(angularSpacing)) - dHole) * 1e3;
                % --- Calculate velocities ---
                vHole = actualDm / (rho * cdHole * totalHoles * holeArea);
                vAnnulus = otherFluidDm / (otherFluidRho * pi * ((rOut)^2 - (rPintle_mm * 1e-3)^2));
                % --- Calculate performance parameters ---
                TMR = (actualDm * vHole) / (otherFluidDm * vAnnulus);
                alpha = rad2deg(acos(1 / (1 + TMR)));
                BF = (holesPerRow * dHole) / (2 * pi * rPintle_mm * 1e-3);
                % --- Store results ---
                results.(configs{configType}){i, j, k} = struct(...
                    'holeDiameter_mm', dHole * 1000, ...
                    'holeArea_mm2', holeArea * 1e6, ...
                    'rowCount', numRows, ...
                    'holesPerRow', holesPerRow, ...
                    'totalHoles', totalHoles, ...
                    'massFlow_kg_s', actualDm, ...
                    'massFlowError_percent', abs(actualDm - dm) / dm * 100, ...
                    'angularSpacing_deg', angularSpacing, ...
                    'arcDistance_mm', arcDistance, ...
                    'rPintle_mm', rPintle_mm, ...
                    'annulusWidth_mm', (rOut * 1e3) - rPintle_mm, ...
                    'vHole_m_s', vHole, ...
                    'vAnnulus_m_s', vAnnulus, ...
                    'TMR', TMR, ...
                    'sprayAngle_deg', alpha, ...
                    'blockageFactor', BF);
            end
        end
    end
end
%% --- Create Analysis Excel File ---
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
analysisFile = sprintf('PintleInjector_Analysis_%s.xlsx', timestamp);
solidworksFile = sprintf('PintleInjector_SolidWorks_%s.xlsx', timestamp);
% --- Input Parameters for Excel ---
inputParams = {'INPUT PARAMETERS', '', ''; ...
    'Parameter', 'Value', 'Units'; ...
    'Injection Pressure (Fuel)', pInj_fuel, 'bar'; ...
    'Injection Pressure (Oxidizer)', pInj_ox, 'bar'; ...
    'Chamber Pressure', pComb, 'bar'; ...
    'Chamber Diameter', dComb, 'mm'; ...
    'O/F Ratio', OF, '-'; ...
    'Oxidizer Mass Flow', dmOx, 'kg/s'; ...
    'Fuel Mass Flow', dmFuel, 'kg/s'; ...
    'Oxidizer Temperature', oxTemp, 'K'; ...
    'Tank Temperature', tankTemperature, 'K'};
inputParams = [inputParams, cell(size(inputParams,1), 13)]; % Pad inputParams to 16 columns
% --- Results Headers for Excel ---
headers = {'Configuration', 'Hole Diameter (mm)', 'Rows', 'Holes/Row', ...
    'Total Holes', 'Pintle Radius (mm)', 'Angular Spacing (deg)', ...
    'Arc Distance (mm)', 'Annulus Width (mm)', 'Mass Flow (kg/s)', ...
    'Flow Error (%)', 'Hole Velocity (m/s)', 'Annulus Velocity (m/s)', ...
    'TMR', 'Spray Angle (deg)', 'Blockage Factor'};
% --- Compile results for Excel in a structured format ---
excelData = {};
excelData = [excelData; inputParams; cell(1,16)]; % Add Input Parameters and a blank