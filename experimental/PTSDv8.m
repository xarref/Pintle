clc; clear

%% Input Parameters
% Operating Conditions
pInj_fuel = 60 * 0.95 - 6;  % Fuel injection pressure [bar]
pInj_ox = 50 * 0.95;        % Oxidizer injection pressure [bar]
pComb = 35;                 % Chamber pressure [bar]
dComb = 76;                 % Chamber diameter [mm]
OF = 2.9;                   % Oxidizer to fuel ratio [-]

% Mass Flow Parameters
dmOx = 1.3502;             % Oxidizer mass flow [kg/s]
dmFuel = dmOx/OF;          % Fuel mass flow [kg/s]

% Fluid Properties
oxidizer = 'N2O';          
fuel = 'Ethanol';
tankTemperature = 298;     % Tank temperature [K]
oxTemp = 273.15;           % Oxidizer temperature [K]

% Discharge Coefficients
cdAnnulus = 0.65;          % Annular injection Cd [-]
cdHole = 0.611;            % Discrete holes Cd [-]

% Geometric Parameters
maxChamberToPintleRatio = 5;
minChamberToPintleRatio = 3;
dPintle = linspace(dComb/maxChamberToPintleRatio, dComb/minChamberToPintleRatio, 5);
rPintle = dPintle/2;

% Configuration Options
dHoleOptions = [0.5, 0.6, 0.8, 1, 1.2] * 1e-3;  % [m]
rowOptions = 1:5;

%% Fluid Properties Calculation
% Get fluid properties using CoolProp
rhoOx = PropSI('D', 'P', pInj_ox*1e5, 'T', oxTemp, oxidizer);
rhoFuel = PropSI('D', 'P', (pInj_fuel-6)*1e5, 'T', 490, fuel);

% Calculate required injection areas
aOx = dmOx/(cdAnnulus*sqrt(2*rhoOx*(pInj_ox-pComb)*1e5));
aFuel = dmFuel/(cdHole*sqrt(2*rhoFuel*(pInj_fuel-pComb)*1e5));

%% Initialize Storage Arrays
configs = {'Fuel_Internal', 'Oxidizer_Internal'};
nConfigs = length(configs);
nHoles = length(dHoleOptions);
nRows = length(rowOptions);
nPintles = length(rPintle);

% Structure to store all results
results = struct();
% Create fields of the structure before the loop.
for c = 1:nConfigs
    results.(configs{c}) = cell(nHoles, nRows, nPintles);
end

%% Calculate Configurations
for configType = 1:nConfigs
    isFuelInternal = (configType == 1);
    
    % Set properties based on configuration
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
        otherFluidRho = rhoFuel;
    end
    
    % Calculate for each combination
    for i = 1:nHoles
        dHole = dHoleOptions(i);
        holeArea = pi*(dHole/2)^2;
        
        for j = 1:nRows
            nRows = rowOptions(j);
            holesPerRow = floor(dm/(cdHole*holeArea*sqrt(2*rho*(pInj-pComb)*1e5)*nRows));
            totalHoles = holesPerRow * nRows;
            actualDm = cdHole * totalHoles * holeArea * sqrt(2*rho*(pInj-pComb)*1e5);
            angularSpacing = 360/holesPerRow;
            
            for k = 1:nPintles
                % Calculate geometrical parameters
                rOut = sqrt(isFuelInternal * aFuel/pi + (~isFuelInternal * aOx/pi) + (rPintle(k)*1e-3)^2);
                arcDistance = ((rPintle(k)*1e-3 * deg2rad(angularSpacing)) - dHole)*1e3;
                
                % Calculate velocities
                vHole = actualDm/(rho * cdHole * totalHoles * holeArea);
                vAnnulus = otherFluidDm/(otherFluidRho * pi * ((rOut)^2 - (rPintle(k)*1e-3)^2));
                
                % Calculate performance parameters
                TMR = (actualDm * vHole)/(otherFluidDm * vAnnulus);
                alpha = rad2deg(acos(1/(1 + TMR)));
                BF = (holesPerRow * dHole)/(2 * pi * rPintle(k)*1e-3);
                
                % Store results
                results.(configs{configType}){i,j,k} = struct(...
                    'holeDiameter', dHole * 1000, ...
                    'holeArea', holeArea * 1e6, ...
                    'rowCount', nRows, ...
                    'holesPerRow', holesPerRow, ...
                    'totalHoles', totalHoles, ...
                    'massFlow', actualDm, ...
                    'massFlowError', abs(actualDm - dm)/dm * 100, ...
                    'angularSpacing', angularSpacing, ...
                    'arcDistance', arcDistance, ...
                    'rPintle', rPintle(k), ...
                    'annulusWidth', (rOut * 1e3) - rPintle(k), ...
                    'vHole', vHole, ...
                    'vAnnulus', vAnnulus, ...
                    'TMR', TMR, ...
                    'alpha', alpha, ...
                    'BF', BF);
            end
        end
    end
end

%% Create Analysis Excel File
analysisFile = 'PintleInjector_Analysis.xlsx';

% Write input parameters
inputParams = {'Parameter', 'Value', 'Units'; ...
    'Injection Pressure (Fuel)', pInj_fuel, 'bar'; ...
    'Injection Pressure (Oxidizer)', pInj_ox, 'bar'; ...
    'Chamber Pressure', pComb, 'bar'; ...
    'Chamber Diameter', dComb, 'mm'; ...
    'O/F Ratio', OF, '-'; ...
    'Oxidizer Mass Flow', dmOx, 'kg/s'; ...
    'Fuel Mass Flow', dmFuel, 'kg/s'; ...
    'Oxidizer Temperature', oxTemp, 'K'; ...
    'Tank Temperature', tankTemperature, 'K'};

writecell({'INPUT PARAMETERS'}, analysisFile, 'Sheet', 'Inputs', 'Range', 'A1');
writecell(inputParams, analysisFile, 'Sheet', 'Inputs', 'Range', 'A3');

% Prepare results headers
headers = {'Hole Diameter (mm)', 'Rows', 'Holes/Row', ...
    'Total Holes', 'Pintle Radius (mm)', 'Angular Spacing (deg)', ...
    'Arc Distance (mm)', 'Annulus Width (mm)', 'Mass Flow (kg/s)', ...
    'Flow Error (%)', 'Hole Velocity (m/s)', 'Annulus Velocity (m/s)', ...
    'TMR', 'Spray Angle (deg)', 'Blockage Factor'};

% Compile all results
allResults = {};
for c = 1:nConfigs
    for i = 1:nHoles
        for j = 1:nRows
            for k = 1:nPintles
                r = results.(configs{c}){i,j,k};
                allResults(end+1,:) = {r.holeDiameter, r.rowCount, ...
                    r.holesPerRow, r.totalHoles, r.rPintle, r.angularSpacing, ...
                    r.arcDistance, r.annulusWidth, r.massFlow, r.massFlowError, ...
                    r.vHole, r.vAnnulus, r.TMR, r.alpha, r.BF};
            end
        end
    end
end

% Write results with parallel comparison
writecell({'COMPREHENSIVE RESULTS'}, analysisFile, 'Sheet', 'Results', 'Range', 'A1');
writecell({'Fuel Internal Configuration'}, analysisFile, 'Sheet', 'Results', 'Range', 'A3');
writecell(headers, analysisFile, 'Sheet', 'Results', 'Range', 'A4');
writecell(allResults(1:size(allResults,1)/2,:), analysisFile, 'Sheet', 'Results', 'Range', 'A5');

writecell({'Oxidizer Internal Configuration'}, analysisFile, 'Sheet', 'Results', 'Range', 'P3');
writecell(headers, analysisFile, 'Sheet', 'Results', 'Range', 'P4');
writecell(allResults(size(allResults,1)/2+1:end,:), analysisFile, 'Sheet', 'Results', 'Range', 'P5');

%% Format Excel File
formatExcelFile(analysisFile);

fprintf('Analysis complete. Files generated:\n');
fprintf('1. %s - Complete analysis\n', analysisFile);

%% Helper Functions
function result = PropSI(output, input1, value1, input2, value2, fluid)
    result = py.CoolProp.CoolProp.PropsSI(output, input1, value1, input2, value2, fluid);
end

function formatExcelFile(filename)
    try
        excel = actxserver('Excel.Application');
        workbook = excel.Workbooks.Open(fullfile(pwd, filename));
        
        % Format each sheet
        for i = 1:workbook.Sheets.Count
            sheet = workbook.Sheets.Item(i);
            sheet.Columns.Item(1).ColumnWidth = 25;
            sheet.Columns.Item(2).ColumnWidth = 15;
            sheet.Columns.Item(3).ColumnWidth = 10;
            if sheet.Columns.Count >= 4
                sheet.Columns.Item(4).ColumnWidth = 40;
            end
            
            % Apply borders and colors
            range = sheet.Range('A1:Z1000');
            range.Borders.Item('xlEdgeBottom').LineStyle = 1;
            range.Borders.Item('xlEdgeTop').LineStyle = 1;
            range.Borders.Item('xlEdgeLeft').LineStyle = 1;
            range.Borders.Item('xlEdgeRight').LineStyle = 1;
            range.Borders.Item('xlInsideVertical').LineStyle = 1;
            range.Borders.Item('xlInsideHorizontal').LineStyle = 1;
            
            % Color headers
            headerRange = sheet.Range('A3:O3');
            headerRange.Interior.Color = rgb2com([173, 216, 230]); % Light blue
            headerRange.Font.Bold = true;
            
            headerRange = sheet.Range('P3:Z3');
            headerRange.Interior.Color = rgb2com([144, 238, 144]); % Light green
            headerRange.Font.Bold = true;
        end
        
        workbook.Save;
        workbook.Close;
        excel.Quit;
        excel.delete;
    catch
        warning('Excel formatting failed. Files saved without formatting.');
    end
end

function comColor = rgb2com(rgbColor)
    % Convert RGB to Excel color format
    comColor = rgbColor(1) + (rgbColor(2) * 256) + (rgbColor(3) * 256 * 256);
end