function pintleInjectorDesignGUI()
    % Create the main figure
    fig = uifigure('Name', 'Pintle Injector Design Tool V5', 'Position', [100 100 1200 800]);
    
    % Input Panel
    inputPanel = uipanel(fig, 'Title', 'Input Parameters', 'Position', [20 650 1160 140]);
    
    % Create input fields with labels
    inputFields = createInputFields(inputPanel);
    
    % Tabbed Results Display
    tabgroup = uitabgroup('Parent', uipanel(fig, 'Title', 'Design Results', 'Position', [20 50 1160 580]), ...
        'Position', [10 10 1140 550]);
    tab1 = uitab(tabgroup, 'Title', 'Geometric Parameters');
    tab2 = uitab(tabgroup, 'Title', 'Injection Characteristics');
    
    % Tables for each tab
    geomTable1 = uitable(tab1, 'Position', [10 10 530 500]);
    geomTable2 = uitable(tab1, 'Position', [550 10 530 500]);
    injectionTable1 = uitable(tab2, 'Position', [10 10 530 500]);
    injectionTable2 = uitable(tab2, 'Position', [550 10 530 500]);
    
    % Compute Button
    uibutton(inputPanel, 'Text', 'Compute', 'Position', [1080 20 80 100], ...
        'ButtonPushedFcn', @(btn,event) computeDesign());
    
    function computeDesign()
        try
            % Collect input values
            inputs = collectInputValues(inputFields);
            
            % Determine configurations to analyze
            configs = determineConfigurations(inputs.configSelection);
            
            % Initialize results storage
            results = cell(numel(configs), 1);
            
            % Run analysis for selected configurations
            for idx = 1:numel(configs)
                configuration = configs(idx);
                
                % Perform calculations
                results{idx} = performPintleDesignAnalysis(...
                    inputs.dmOx, inputs.OF, inputs.pInj, inputs.pComb, ...
                    inputs.dComb, inputs.chamberPintleRatio, inputs.pintleAngle, configuration);
            end
            
            % Update tables with results
            updateResultsTables(results, geomTable1, geomTable2, injectionTable1, injectionTable2);
            
        catch ME
            % Display any errors
            uialert(fig, getReport(ME), 'Calculation Error');
        end
    end
end

function inputFields = createInputFields(inputPanel)
    % Mass Flow Rate (Oxidizer)
    inputFields.dmOxLabel = uilabel(inputPanel, 'Text', 'Oxidizer Mass Flow Rate (kg/s)', 'Position', [10 90 200 22]);
    inputFields.dmOxEdit = uieditfield(inputPanel, 'numeric', 'Position', [220 90 100 22], 'Value', 1.5);
    
    % OF Ratio
    inputFields.ofLabel = uilabel(inputPanel, 'Text', 'Oxidizer/Fuel Ratio', 'Position', [10 60 200 22]);
    inputFields.ofEdit = uieditfield(inputPanel, 'numeric', 'Position', [220 60 100 22], 'Value', 3.8);
    
    % Injection Pressures
    inputFields.pInjLabel = uilabel(inputPanel, 'Text', 'Injection Pressure (bar)', 'Position', [350 90 200 22]);
    inputFields.pInjEdit = uieditfield(inputPanel, 'numeric', 'Position', [550 90 100 22], 'Value', 51);
    
    % Combustion Chamber Pressure
    inputFields.pCombLabel = uilabel(inputPanel, 'Text', 'Combustion Chamber Pressure (bar)', 'Position', [350 60 200 22]);
    inputFields.pCombEdit = uieditfield(inputPanel, 'numeric', 'Position', [550 60 100 22], 'Value', 35);
    
    % Chamber Diameter
    inputFields.dCombLabel = uilabel(inputPanel, 'Text', 'Chamber Diameter (mm)', 'Position', [350 30 200 22]);
    inputFields.dCombEdit = uieditfield(inputPanel, 'numeric', 'Position', [550 30 100 22], 'Value', 100);
    
    % Specific Chamber/Pintle Ratio
    inputFields.chamberPintleRatioLabel = uilabel(inputPanel, 'Text', 'Specific Chamber/Pintle Ratio', 'Position', [10 30 200 22]);
    inputFields.chamberPintleRatioEdit = uieditfield(inputPanel, 'numeric', 'Position', [220 30 100 22], 'Value', 4);
    
    % Pintle Angle
    inputFields.pintleAngleLabel = uilabel(inputPanel, 'Text', 'Pintle Angle (degrees)', 'Position', [700 30 200 22]);
    inputFields.pintleAngleEdit = uieditfield(inputPanel, 'numeric', 'Position', [900 30 100 22], 'Value', 30);
    
    % Configuration Selection
    inputFields.configLabel = uilabel(inputPanel, 'Text', 'Configuration', 'Position', [700 60 200 22]);
    inputFields.configDropdown = uidropdown(inputPanel, 'Items', {'Both', 'Ox Internal, Fuel External', 'Fuel Central, Ox External'}, ...
        'Position', [900 60 200 22], 'Value', 'Both');
end

function inputs = collectInputValues(inputFields)
    inputs.dmOx = inputFields.dmOxEdit.Value;
    inputs.OF = inputFields.ofEdit.Value;
    inputs.pInj = inputFields.pInjEdit.Value;
    inputs.pComb = inputFields.pCombEdit.Value;
    inputs.dComb = inputFields.dCombEdit.Value;
    inputs.chamberPintleRatio = inputFields.chamberPintleRatioEdit.Value;
    inputs.pintleAngle = inputFields.pintleAngleEdit.Value;
    inputs.configSelection = inputFields.configDropdown.Value;
end

function configs = determineConfigurations(configSelection)
    switch configSelection
        case 'Both'
            configs = [1, 2];
        case 'Ox Internal, Fuel External'
            configs = 1;
        case 'Fuel Central, Ox External'
            configs = 2;
    end
end

function updateResultsTables(results, geomTable1, geomTable2, injectionTable1, injectionTable2)
% Prepare data for tables
geomDataFirst = cell(0,2);
geomDataSecond = cell(0,2);
injectionDataFirst = cell(0,2);
injectionDataSecond = cell(0,2);

for idx = 1:numel(results)
    res = results{idx};

    if res.configuration == 1
        configName = 'Ox Internal, Fuel External';
    else
        configName = 'Fuel Central, Ox External';
    end

    % Geometric Parameters
    geomDataCurrent = {...
        'Configuration', configName;
        'Minimum Injection Area (aMin)', sprintf('%.4f mm²', res.aMin);
        'Fuel Injection Area', sprintf('%.4f mm²', res.aFuel);
        'Cavity Area (aCg)', sprintf('%.4f mm²', res.aCg);
        'Pintle Outer Radius (rPost)', sprintf('%.4f mm', res.rPost);
        'Pintle Support Radius (rPr)', sprintf('%.4f mm', res.rPr);
        'Second Channel Radius', sprintf('%.4f mm', res.rSecondChannel);
        'First Channel Wall Thickness', sprintf('%.4f mm', res.twall_first_channel);
        'External Annulus Area', sprintf('%12.4f mm²', res.externalAnnulusArea);
        'External Annulus Width', sprintf('%.4f mm', res.externalAnnulusWidth);
        'Channel Open Length (LOpen)', sprintf('%.4f mm', res.LOpen); % Added parameter
        'Minimum Length (LMin)', sprintf('%.4f mm', res.LMin) % Added parameter
        };

    % Injection Characteristics
    injectionDataCurrent = {...
        'Configuration', configName;
        'Oxidizer Mass Flow Rate', sprintf('%.4f kg/s', res.dmOx);
        'Fuel Mass Flow Rate', sprintf('%.4f kg/s', res.dmFuel);
        'Total Mass Flow Rate', sprintf('%.4f kg/s', res.dmOx + res.dmFuel);
        'Injection Velocity (Oxidizer)', sprintf('%.4f m/s', res.oxInjectionVelocity);
        'Injection Velocity (Fuel)', sprintf('%.4f m/s', res.fuelInjectionVelocity);
        'Droplet Sauter Mean Diameter', sprintf('%.4f μm', res.smdDropletSize)
        };

    % Assign to appropriate table
    if idx == 1
        geomDataFirst = geomDataCurrent;
        injectionDataFirst = injectionDataCurrent;
    else
        geomDataSecond = geomDataCurrent;
        injectionDataSecond = injectionDataCurrent;
    end

end

% Set table data
set(geomTable1, 'Data', geomDataFirst, 'ColumnName', {'Parameter', 'Value'});
set(geomTable2, 'Data', geomDataSecond, 'ColumnName', {'Parameter', 'Value'});
set(injectionTable1, 'Data', injectionDataFirst, 'ColumnName', {'Parameter', 'Value'});
set(injectionTable2, 'Data', injectionDataSecond, 'ColumnName', {'Parameter', 'Value'});

end
function result = performPintleDesignAnalysis(dmOx, OF, pInj, pComb, dComb, chamberPintleRatio, pintleAngle, configuration)
    % Oxidizer and Fuel Properties
    oxidizer = 'N2O';
    fuel = 'Ethanol';
    
    % Tank temperature
    tankTemperature = 298;
    
    % Discharge coefficients
    cdOx = 0.65;
    cdFuelHole = 0.6;
    
    % Mass flow rates
    dmFuel = dmOx/OF;
    
    % Compute fuel properties
    rhoFuelPropSI = PropSI('D', 'P', pInj*1e5*0.8, 'T', 450, fuel); %density of ethanol at 450K, assuming 20% pressure drop in channels (temperature is within the safe margin of the specific heat of ethanol)
    
    % Compute pintle diameter range
    dPintle = dComb/chamberPintleRatio;
  
    
    % NHNE MODEL calculations
    % Interpolating liquid and vapor properties prior injection
    hLi = PropSI('H', 'P', pInj*1e5, 'Q', 0, oxidizer);
    rhoLi = PropSI('D', 'P', pInj*1e5, 'Q', 0, oxidizer);
    sL = PropSI('S', 'P', pInj*1e5, 'Q', 0, oxidizer);
    
    % Load choke parameter function
    fChockParamFun = load("fChockParamFun.mat");
    poOx = fChockParamFun.f_chock_param(pInj) * pInj;
    
    % More interpolations
    hLo = PropSI('H', 'P', poOx*1e5, 'Q', 0, oxidizer);
    hVo = PropSI('H', 'P', poOx*1e5, 'Q', 1, oxidizer);
    sLo = PropSI('S', 'P', poOx*1e5, 'Q', 0, oxidizer);
    sVo = PropSI('S', 'P', poOx*1e5, 'Q', 1, oxidizer);
    
    rhoLo = PropSI('D', 'P', poOx*1e5, 'Q', 0, oxidizer);
    rhoVo = PropSI('D', 'P', poOx*1e5, 'Q', 1, oxidizer);
    
    % Compute quality and other parameters
    so = sL;
    Xo = (so - sLo)/(sVo - sLo);
    ho = hLo * (1 - Xo) + hVo * Xo;
    rhoo = rhoVo * rhoLo / (rhoLo * Xo + rhoVo * (1 - Xo));
    
    % Compute vapor pressure
    PV = PropSI('P', 'T', tankTemperature, 'Q', 1, oxidizer) / 1e5;
    
    % Compute k parameter
    k = sqrt((pInj - poOx)/(PV - poOx));
    
    % Compute injection areas using NHNE and SPI models
    aOx = dmOx/( (1-(1/(1+k)))*cdOx*sqrt(2*rhoLi*(pInj-poOx)*1e5) ...
        + (1/(1+k))*cdOx*rhoo*sqrt(2*(hLi-ho)) ); %in m 
    aOx = aOx * 1e6; % Convert to mm²
   
    
    % Fuel injection area
    aFuel = dmFuel/(cdFuelHole*sqrt(2*rhoFuelPropSI*(pInj-pComb)*1e5));
    aFuel = aFuel * 1e6; % Convert to mm²
    
    % Determine minimum area based on configuration
    if configuration == 1
        aMin = aOx;  % Ox Internal, Fuel External
        externalAnnulusArea = aFuel;
    else
        aMin = aFuel; % Fuel Central, Ox External
        externalAnnulusArea = aOx;
    end
    
    % Design calculations
    rPr_fixed = 4; % Minimum pintle root radius
    ratio_AcgAmin = 1.5;
    aCg = ratio_AcgAmin * aMin;
    
    rPost = dPintle/2;
    rPr = rPr_fixed;
    rCg = sqrt((aCg/pi) + rPr^2);
    twall_first_channel = rPost - rCg;
    
    % Compute second channel radius (configuration-dependent)
    if configuration == 1 % Ox Internal, Fuel External
        rSecondChannel = sqrt((aFuel/pi) + rPost^2);
    else % Fuel Central, Ox External
        rSecondChannel = sqrt((aOx/pi) + rPost^2);
    end
    
    % Convert pintle angle to radians
    pintleAngleRad = deg2rad(pintleAngle);
    
    % External Annulus Area Calculation
    externalAnnulusWidth = rPost * tan(pintleAngleRad/2);
    
    
    % Compute channel open length and minimum length
    % Note: This is a simplified version. In the full code, multiple angles are considered
    % For this implementation, we'll use the provided pintle angle
    LOpen = ( (rPost-1.5) - ( sqrt((rPost - 1.5)^2 - aMin*(sin(pintleAngleRad)/pi) ) ) )/ (sin(pintleAngleRad)*cos(pintleAngleRad));
    LMin = LOpen * cos(pintleAngleRad);
    
    % Simple injection velocity estimation
    oxInjectionVelocity = sqrt(2 * (pInj - poOx) * 1e5 / rhoLi);
    fuelInjectionVelocity = sqrt(2 * (pInj - pComb) * 1e5 / rhoFuelPropSI);
    
    % Estimated Sauter Mean Diameter (SMD) using a simplified correlation
    We = (rhoLi * oxInjectionVelocity^2 * externalAnnulusWidth) / (0.072); % Weber number
    Oh = (sqrt(0.072 * rhoLi) * (1e-3)^0.5) / (sqrt(0.072/rhoLi)); % Ohnesorge number
    smdDropletSize = 0.5 * externalAnnulusWidth * (We^0.5 / Oh^0.5);
    
    % Package and return results
    result = struct(...
        'configuration', configuration, ...
        'aMin', aMin, ...
        'aFuel', (configuration == 1) * aFuel + (configuration ~= 1) * aOx, ...
        'aCg', aCg, ...
        'rPost', rPost, ...
        'rPr', rPr, ...
        'rSecondChannel', rSecondChannel, ...
        'twall_first_channel', twall_first_channel, ...
        'externalAnnulusArea', externalAnnulusArea, ...
        'externalAnnulusWidth', externalAnnulusWidth, ...
        'LOpen', LOpen, ...
        'LMin', LMin, ...
        'dmOx', dmOx, ...
        'dmFuel', dmFuel, ...
        'oxInjectionVelocity', oxInjectionVelocity, ...
        'fuelInjectionVelocity', fuelInjectionVelocity, ...
        'smdDropletSize', smdDropletSize ...
    );
end

function result = PropSI(output, input1, value1, input2, value2, fluid)
    % Wrapper function for CoolProp PropsSI
    result = py.CoolProp.CoolProp.PropsSI(output, input1, value1, input2, value2, fluid);
end

% Main function to run the GUI
function main()
    pintleInjectorDesignGUI();
end

% Execute the GUI
main();