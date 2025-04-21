clc; clear

%% Input Parameters
% Operating Conditions
pInj_fuel = 60 * 0.95 - 6;  % Fuel injection pressure [bar]
pInj_ox = 50 * 0.95;        % Oxidizer injection pressure [bar]
pComb = 35;                 % Chamber pressure [bar]
dComb = 90.23;                 % Chamber diameter [mm]
%dComb = 76;                 % Chamber diameter [mm]
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
                rOut = sqrt(isFuelInternal * aOx/pi + (~isFuelInternal * aFuel/pi) + (rPintle(k)*1e-3)^2);
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

%% Plots design parameters

for c = 1:nConfigs
    if c == 1
        fprintf('-----------------------DESIGN PARAMETERS-------------------------------\n');
        fprintf('-----------------------------------------------------------\n');
        fprintf('Fuel Internal with Holes and N2O Annulus \n\n\n');
        fprintf('Computed fuel injection area: %.4f mm²\n', aFuel*1e6);
        fprintf('Computer oxidizer injection area: %.4f mm²\n', aOx*1e6)
        fprintf('-----------------------------------------------------------\n');

        fprintf('Fuel Discrete Hole Configuration Analysis\n\n');
        fprintf('------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');
        fprintf('| Hole Dia. | Hole Area  | Rows | Holes/Row   | Total Holes  | Mass Flow (kg/s) | Error (%%) || Angular Delta | Arc Distance (mm)  | rPintle (mm) | Annulus Width (mm)  |\n');
        fprintf('|-----------------------------------------------------------------------------------------------------------------------------------------------------------------|\n');
        for i = 1:nHoles
            for j = 1:nRows
                for k = 1:nPintles

                    fprintf('| %.2f mm   | %.4f mm² |   %d  |     %4d    |     %4d     |      %.4f      |   %2.2f%%   ||   %6.2f deg  |     %6.3f mm      |   %6.2f mm   | %6.2f mm \n', ...
                        results.Fuel_Internal{i,j,k}.holeDiameter, ...
                        results.Fuel_Internal{i,j,k}.holeArea, ...
                        results.Fuel_Internal{i,j,k}.rowCount, ...
                        results.Fuel_Internal{i,j,k}.holesPerRow, ...
                        results.Fuel_Internal{i,j,k}.totalHoles, ...
                        results.Fuel_Internal{i,j,k}.massFlow, ...
                        results.Fuel_Internal{i,j,k}.massFlowError, ...
                        results.Fuel_Internal{i,j,k}.angularSpacing, ...
                        results.Fuel_Internal{i,j,k}.arcDistance, ...
                        results.Fuel_Internal{i,j,k}.rPintle, ...
                        results.Fuel_Internal{i,j,k}.annulusWidth ...
                        );
                end
            end
             fprintf('|------------------------------------------------------------------------------------------------------------------------------------------------------------------|\n\n\n\n\n');

        end

      else

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

          for i = 1:nHoles
            for j = 1:nRows
                for k = 1:nPintles

                    fprintf('| %.2f mm   | %.4f mm² |   %d  |     %4d    |     %4d     |      %.4f      |   %2.2f%%   ||   %6.2f deg  |     %6.3f mm      |   %6.2f mm   | %6.2f mm \n', ...
                        results.Oxidizer_Internal{i,j,k}.holeDiameter, ...
                        results.Oxidizer_Internal{i,j,k}.holeArea, ...
                        results.Oxidizer_Internal{i,j,k}.rowCount, ...
                        results.Oxidizer_Internal{i,j,k}.holesPerRow, ...
                        results.Oxidizer_Internal{i,j,k}.totalHoles, ...
                        results.Oxidizer_Internal{i,j,k}.massFlow, ...
                        results.Oxidizer_Internal{i,j,k}.massFlowError, ...
                        results.Oxidizer_Internal{i,j,k}.angularSpacing, ...
                        results.Oxidizer_Internal{i,j,k}.arcDistance, ...
                        results.Oxidizer_Internal{i,j,k}.rPintle, ...
                        results.Oxidizer_Internal{i,j,k}.annulusWidth ...
                        );
                end
            end
             fprintf('|------------------------------------------------------------------------------------------------------------------------------------------------------------------|\n\n\n\n\n');

          end
    end
end

%% Plot Performances

for c = 1:nConfigs
    if c == 1
        fprintf('Fuel Internal with Holes and N2O Annulus \n\n\n');
        fprintf('--------------------------------------------------PERFORMANCE PARAMETERS--------------------------------------------------\n');
        fprintf('| Hole Dia. | Hole Area  | rPintle (mm)  |   v_oxid (m/s)  |   v_fuel (m/s) |    TMR     |   alpha   |     BF     | Rows |  Holes/Row  | Total Holes |  \n');
         fprintf('|--------------------------------------------------------------------------------------------------------------------------------------------------------|\n');

        for i = 1:nHoles
            for j = 1:nRows
                for k = 1:nPintles

                    fprintf('| %.2f mm   | %.4f mm² |   %6.2f mm   |   %7.4f m/s   |   %7.4f m/s  |   %6.4f   |   %5.2f°  |   %6.4f   |   %d  |    %4d     |    %4d    |   \n', ...
                        results.Fuel_Internal{i,j,k}.holeDiameter, ...
                        results.Fuel_Internal{i,j,k}.holeArea, ...
                        results.Fuel_Internal{i,j,k}.rPintle, ...
                        results.Fuel_Internal{i,j,k}.vAnnulus, ...
                        results.Fuel_Internal{i,j,k}.vHole, ...
                        results.Fuel_Internal{i,j,k}.TMR, ...
                        results.Fuel_Internal{i,j,k}.alpha, ...
                        results.Fuel_Internal{i,j,k}.BF, ...
                        results.Fuel_Internal{i,j,k}.rowCount, ...
                        results.Fuel_Internal{i,j,k}.holesPerRow, ...
                        results.Fuel_Internal{i,j,k}.totalHoles ...
                        );
                end
            end
             fprintf('|------------------------------------------------------------------------------------------------------------------------------------------------------------------|\n\n\n\n\n');

        end

      else
          fprintf('N2O Internal with Holes and fuel Annulus \n\n\n');
          fprintf('--------------------------------------------------PERFORMANCE PARAMETERS--------------------------------------------------\n');
          fprintf('| Hole Dia. | Hole Area  | rPintle (mm)  |   v_oxid (m/s)  |   v_fuel (m/s) |    TMR     |   alpha   |     BF     | Rows |  Holes/Row  | Total Holes |  \n');
           fprintf('|--------------------------------------------------------------------------------------------------------------------------------------------------------|\n');

          for i = 1:nHoles
            for j = 1:nRows
                for k = 1:nPintles

                    fprintf('| %.2f mm   | %.4f mm² |   %6.2f mm   |   %7.4f m/s   |   %7.4f m/s  |   %6.4f   |   %5.2f°  |   %6.4f   |   %d  |    %4d     |    %4d    |   \n', ...
                        results.Oxidizer_Internal{i,j,k}.holeDiameter, ...
                        results.Oxidizer_Internal{i,j,k}.holeArea, ...
                        results.Oxidizer_Internal{i,j,k}.rPintle, ...
                        results.Oxidizer_Internal{i,j,k}.vHole, ...
                        results.Oxidizer_Internal{i,j,k}.vAnnulus, ...
                        results.Oxidizer_Internal{i,j,k}.TMR, ...
                        results.Oxidizer_Internal{i,j,k}.alpha, ...
                        results.Oxidizer_Internal{i,j,k}.BF, ...
                        results.Oxidizer_Internal{i,j,k}.rowCount, ...
                        results.Oxidizer_Internal{i,j,k}.holesPerRow, ...
                        results.Oxidizer_Internal{i,j,k}.totalHoles ...
                        );
                end
            end
             fprintf('|------------------------------------------------------------------------------------------------------------------------------------------------------------------|\n\n\n\n\n');

          end
    end
end

%% Add Timestamp to File Names
timestamp = datestr(now, 'yyyy-mm-dd');
analysisFile = sprintf('%s_PintleInjector_Analysis.xlsx', timestamp);
solidworksFile = sprintf('%s_PintleInjector_SolidWorks.xlsx', timestamp);

%% Create Analysis Excel File
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
headers = {'Configuration', 'Hole Diameter (mm)', 'Rows', 'Pintle Radius (mm)', ...
    'Holes/Row', 'Total Holes', 'Angular Spacing (deg)', ...
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
                allResults(end+1,:) = {configs{c}, r.holeDiameter, r.rowCount, ...
                    r.rPintle, r.holesPerRow, r.totalHoles, r.angularSpacing, ...
                    r.arcDistance, r.annulusWidth, r.massFlow, r.massFlowError, ...
                    r.vHole, r.vAnnulus, r.TMR, r.alpha, r.BF};
            end
        end
    end
end

% Write results
writecell({'COMPREHENSIVE RESULTS'}, analysisFile, 'Sheet', 'Results', 'Range', 'A1');
writecell(headers, analysisFile, 'Sheet', 'Results', 'Range', 'A3');
writecell(allResults, analysisFile, 'Sheet', 'Results', 'Range', 'A4');

%% Apply Formatting to Analysis Excel File
try
    excel = actxserver('Excel.Application');
    workbook = excel.Workbooks.Open(fullfile(pwd, analysisFile));
    sheet = workbook.Sheets.Item('Results');
    
    % Enable filters
    sheet.Range('A3:P3').AutoFilter(1);
    
    % Apply background colors and conditional formatting
    lastRow = size(allResults, 1) + 3; % Header rows are 1-3
    for row = 4:lastRow
        config = sheet.Range(sprintf('A%d', row)).Value;
        
        % Set background color based on configuration
        if strcmp(config, 'Fuel_Internal')
            sheet.Rows.Item(row).Interior.Color = rgb2com([255, 182, 193]); % Light red
        elseif strcmp(config, 'Oxidizer_Internal')
            sheet.Rows.Item(row).Interior.Color = rgb2com([173, 216, 230]); % Light blue
        end
        
        % Check for negative values and override background color
        for col = 2:16 % Columns B to P
            cellValue = sheet.Range(sprintf('%s%d', char(64 + col), row)).Value;
            if isnumeric(cellValue) && cellValue < 0
                sheet.Rows.Item(row).Interior.Color = rgb2com([255, 0, 0]); % Dark red
                break; % Highlight the entire row if any value is negative
            end
        end
    end
    
    workbook.Save;
    workbook.Close;
    excel.Quit;
    excel.delete;
catch
    warning('Excel formatting failed. Files saved without formatting.');
end

%% Create SolidWorks Excel File
% Selected configuration (example: 2 rows, 0.6mm holes, Fuel Internal)
selectedConfig = 'Fuel_Internal';
selectedHoleDia = 0.6;
selectedRows = 2;
selectedPintle = 3;  % Middle pintle size

% Find indices
holeDiaIndex = find(dHoleOptions*1000 == selectedHoleDia);
rowIndex = find(rowOptions == selectedRows);
r = results.(selectedConfig){holeDiaIndex,rowIndex,selectedPintle};

% Create SolidWorks parameters file
solidworksParams = {'Parameter', 'Value', 'Units', 'Description'; ...
    % Design Parameters
    'DESIGN_Configuration', selectedConfig, '-', 'Injector configuration type'; ...
    'DESIGN_PintleDiameter', r.rPintle*2, 'mm', 'Pintle external diameter'; ...
    'DESIGN_HoleDiameter', r.holeDiameter, 'mm', 'Injection hole diameter'; ...
    'DESIGN_RowCount', r.rowCount, '-', 'Number of hole rows'; ...
    'DESIGN_HolesPerRow', r.holesPerRow, '-', 'Holes in each row'; ...
    'DESIGN_AngularSpacing', r.angularSpacing, 'deg', 'Angle between holes'; ...
    'DESIGN_RowSpacing', 2.5, 'mm', 'Axial distance between rows'; ...
    'DESIGN_AnnulusWidth', r.annulusWidth, 'mm', 'Width of annular gap'; ...
    
    % Material Properties (Inconel 718)
    'MAT_Density', 8.19, 'g/cm³', 'Density'; ...
    'MAT_TensileStrength', 1375, 'MPa', 'Ultimate tensile strength'; ...
    'MAT_YieldStrength', 1100, 'MPa', 'Yield strength'; ...
    'MAT_ThermalConductivity', 11.4, 'W/m·K', 'Thermal conductivity'; ...
    'MAT_ThermalExpansion', 13, 'μm/m·K', 'Thermal expansion coefficient'; ...
    
    % Manufacturing Constraints (SLM)
    'MFG_MinWallThickness', 0.4, 'mm', 'Minimum wall thickness'; ...
    'MFG_MinHoleDiameter', 0.4, 'mm', 'Minimum hole diameter'; ...
    'MFG_MinFeatureSize', 0.4, 'mm', 'Minimum feature size'; ...
    'MFG_MinRadiusInternal', 0.2, 'mm', 'Minimum internal radius'; ...
    'MFG_MinRadiusExternal', 0.2, 'mm', 'Minimum external radius'; ...
    'MFG_MaxOverhangAngle', 46, 'deg', 'Maximum overhang angle'; ...
    'MFG_SurfaceRoughness', 6.3, 'μm Ra', 'Expected surface finish'; ...
    'MFG_PositionTolerance', 0.1, 'mm', 'Feature position tolerance'; ...
    'MFG_DiameterTolerance', 0.05, 'mm', 'Hole diameter tolerance'};

writecell(solidworksParams, solidworksFile);

%% Format Excel Files
formatExcelFile(analysisFile);
formatExcelFile(solidworksFile);

fprintf('Analysis complete. Files generated:\n');
fprintf('1. %s - Complete analysis\n', analysisFile);
fprintf('2. %s - SolidWorks parameters\n', solidworksFile);

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