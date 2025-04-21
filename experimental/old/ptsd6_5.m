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
analysisFile = 'PintleInjector_Analysis.xlsx';

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
excelData = [excelData; inputParams; cell(1,16)]; % Add Input Parameters and a blank row

for i = 1:nHoles
    for j = 1:nRows
        for k = 1:nPintles
            % --- Parameter combination header ---
            paramHeader = {sprintf('Hole Dia: %.2f mm, Rows: %d, Pintle Radius: %.1f mm', ...
                dHoleOptions(i) * 1000, rowOptions(j), rPintle(k)), '', '', '', '', '', '', '', '', '', '', '', '', '', '', ''};
            excelData = [excelData; paramHeader; headers];

            % --- Results for Fuel Internal and Oxidizer Internal configurations ---
            rFuelInt = results.Fuel_Internal{i, j, k};
            rOxInt = results.Oxidizer_Internal{i, j, k};

            resultRows = [
                {'Fuel Internal', rFuelInt.holeDiameter_mm, rFuelInt.rowCount, rFuelInt.holesPerRow, rFuelInt.totalHoles, rFuelInt.rPintle_mm, rFuelInt.angularSpacing_deg, rFuelInt.arcDistance_mm, rFuelInt.annulusWidth_mm, rFuelInt.massFlow_kg_s, rFuelInt.massFlowError_percent, rFuelInt.vHole_m_s, rFuelInt.vAnnulus_m_s, rFuelInt.TMR, rFuelInt.sprayAngle_deg, rFuelInt.blockageFactor};
                {'Oxidizer Internal', rOxInt.holeDiameter_mm, rOxInt.rowCount, rOxInt.holesPerRow, rOxInt.totalHoles, rOxInt.rPintle_mm, rOxInt.angularSpacing_deg, rOxInt.arcDistance_mm, rOxInt.annulusWidth_mm, rOxInt.massFlow_kg_s, rOxInt.massFlowError_percent, rOxInt.vHole_m_s, rOxInt.vAnnulus_m_s, rOxInt.TMR, rOxInt.sprayAngle_deg, rOxInt.blockageFactor}
                ];
            excelData = [excelData; resultRows; cell(1,16)]; % Add blank row after each parameter set
        end
    end
end

% --- Write data to Excel ---
try
    writecell({'Pintle Injector Analysis Results'}, analysisFile, 'Sheet', 'Results', 'Range', 'A1');
    writecell(excelData, analysisFile, 'Sheet', 'Results', 'Range', 'A2'); % Start writing results below the title.
    fprintf('Analysis results written to: %s\n', analysisFile);
catch
    warning('Failed to write analysis results to Excel file: %s', analysisFile);
end


%% --- Create SolidWorks Excel File ---
solidworksFile = 'PintleInjector_SolidWorks.xlsx';

% --- Selected configuration for SolidWorks parameters (example) ---
selectedConfig = 'Fuel_Internal';
selectedHoleDia = 0.6;
selectedRows = 2;
selectedPintle = 3; % Middle pintle size

holeDiaIndex = find(dHoleOptions * 1000 == selectedHoleDia);
rowIndex = find(rowOptions == selectedRows);


if ~isempty(holeDiaIndex) && ~isempty(rowIndex) && selectedPintle <= length(rPintle)
    r = results.(selectedConfig){holeDiaIndex, rowIndex, selectedPintle};


    solidworksParams = {'Parameter', 'Value', 'Units', 'Description'; ...
        % Design Parameters
        'DESIGN_Configuration', selectedConfig, '-', 'Injector configuration type'; ...
        'DESIGN_PintleDiameter', r.rPintle_mm * 2, 'mm', 'Pintle external diameter'; ...
        'DESIGN_HoleDiameter', r.holeDiameter_mm, 'mm', 'Injection hole diameter'; ...
        'DESIGN_RowCount', r.rowCount, '-', 'Number of hole rows'; ...
        'DESIGN_HolesPerRow', r.holesPerRow, '-', 'Holes in each row'; ...
        'DESIGN_AngularSpacing', r.angularSpacing_deg, 'deg', 'Angle between holes'; ...
        'DESIGN_RowSpacing', 2.5, 'mm', 'Axial distance between rows'; ...
        'DESIGN_AnnulusWidth', r.annulusWidth_mm, 'mm', 'Width of annular gap'; ...

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

    try
        writecell(solidworksParams, solidworksFile);
        fprintf('SolidWorks parameters written to: %s\n', solidworksFile);
    catch
        warning('Failed to write SolidWorks parameters to Excel file: %s', solidworksFile);
    end
else
    warning('Selected configuration parameters for SolidWorks file are invalid. SolidWorks file not generated.');
end


%% --- Format Excel Files ---
% try
%     formatExcelFile(analysisFile);
%     fprintf('Formatted Excel file: %s\n', analysisFile);
% catch
%     warning('Excel formatting failed for: %s', analysisFile);
% end
%
% try
%     formatExcelFile(solidworksFile);
%     fprintf('Formatted Excel file: %s\n', solidworksFile);
% catch
%     warning('Excel formatting failed for: %s', solidworksFile);
% end


fprintf('Analysis complete. Files generated:\n');
fprintf('1. %s - Complete analysis\n', analysisFile);
fprintf('2. %s - SolidWorks parameters (if configuration valid)\n', solidworksFile);


%% --- Helper Functions ---
function result = PropSI(output, input1, value1, input2, value2, fluid)
    % Wrapper for CoolProp PropsSI function with error handling
    try
        result = py.CoolProp.CoolProp.PropsSI(output, input1, value1, input2, value2, fluid);
    catch
        error('CoolProp PropSI function call failed for fluid: %s, output: %s, input1: %s, input2: %s.', fluid, output, input1, input2);
    end
end

% Removed formatExcelFile function
% function formatExcelFile(filename)
%     % Formats the Excel file to improve readability
%     excelApp = [];
%     excelWorkbook = [];
%     try
%         excelApp = actxserver('Excel.Application');
%         excelWorkbook = excelApp.Workbooks.Open(fullfile(pwd, filename), 0, false); % Open workbook, ReadOnly=false (default), Editable=true
%         excelWorkbook.Saved = true; % Prevent "Save As" prompt
%
%         % --- Format Sheets ---
%         for sheetIndex = 1:excelWorkbook.Sheets.Count
%             excelSheet = excelWorkbook.Sheets.Item(sheetIndex);
%
%             % --- Adjust column widths for readability ---
%             excelSheet.Columns.Item(1).ColumnWidth = 25;  % Parameter column
%             excelSheet.Columns.Item(2).ColumnWidth = 18;  % Value column
%             excelSheet.Columns.Item(3).ColumnWidth = 12;  % Units column
%             if excelSheet.Columns.Count >= 4
%                 excelSheet.Columns.Item(4).ColumnWidth = 30; % Description column (if exists)
%             end
%
%             % --- Bold Header Rows ---
%             headerRange = excelSheet.Range('A1:P1');
%             headerRange.Font.Bold = true;
%             headerRange2 = excelSheet.Range('A2:P2');
%             headerRange2.Font.Bold = true;
%
%             % --- Bold Parameter Group Headers and add bottom border ---
%             paramHeaderRows = excelSheet.Range('A3:A1000');
%             paramHeaders = paramHeaderRows.Value;
%             for rowIdx = 1:size(paramHeaders, 1)
%                 if ~isempty(paramHeaders{rowIdx, 1}) && contains(paramHeaders{rowIdx, 1}, 'Hole Dia')
%                     currentRow = excelSheet.Rows.Item(rowIdx + 2);
%                     currentRow.Font.Bold = true;
%
%                     borderRange = excelSheet.Range(['A', num2str(rowIdx + 2), ':', 'P', num2str(rowIdx + 2)]);
%                     borderRange.Borders.Item(Excel.XlBordersIndex.xlEdgeBottom).LineStyle = Excel.XlLineStyle.xlContinuous;
%                     borderRange.Borders.Item(Excel.XlBordersIndex.xlEdgeBottom).Weight = Excel.XlBorderWeight.xlThin;
%                 end
%             end
%         end
%
%         excelWorkbook.Save();
%         excelWorkbook.Close(false); %Modified to prevent save prompt again - should close without saving again as already saved
%         excelApp.Quit();
%         excelApp.delete();
%
%     catch excelError
%         warning('Excel formatting failed for file: %s. Error: %s', filename, excelError.message);
%         if ~isempty(excelWorkbook) && isa(excelWorkbook, 'COM.Excel_TLB.Workbook')
%             excelWorkbook.Close(false);
%         end
%         if ~isempty(excelApp) && isa(excelApp, 'COM.Excel_TLB.Application')
%             excelApp.Quit();
%             excelApp.delete();
%         end
%     end
% end