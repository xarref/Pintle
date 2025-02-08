
%% changelog:
% prima versione della GUI, associato a pintleMasterV2
% non c'è un granchè di roba. 
% todo: sostituire tutto con propSi



function pintleInjectorDesignGUI()
    % Create the main figure
    fig = uifigure('Name', 'Pintle Injector Design Tool', 'Position', [100 100 1000 700]);
    
    % Input Panel
    inputPanel = uipanel(fig, 'Title', 'Input Parameters', 'Position', [20 550 960 140]);
    
    % Mass Flow Rate (Oxidizer)
    uilabel(inputPanel, 'Text', 'Oxidizer Mass Flow Rate (kg/s)', 'Position', [10 90 200 22]);
    dmOxEdit = uieditfield(inputPanel, 'numeric', 'Position', [220 90 100 22], 'Value', 1.3);
    
    % OF Ratio
    uilabel(inputPanel, 'Text', 'Oxidizer/Fuel Ratio', 'Position', [10 60 200 22]);
    ofEdit = uieditfield(inputPanel, 'numeric', 'Position', [220 60 100 22], 'Value', 3.8);
    
    % Injection Pressures
    uilabel(inputPanel, 'Text', 'Injection Pressure (bar)', 'Position', [350 90 200 22]);
    pInjEdit = uieditfield(inputPanel, 'numeric', 'Position', [550 90 100 22], 'Value', 60*0.85);
    
    % Combustion Chamber Pressure
    uilabel(inputPanel, 'Text', 'Combustion Chamber Pressure (bar)', 'Position', [350 60 200 22]);
    pCombEdit = uieditfield(inputPanel, 'numeric', 'Position', [550 60 100 22], 'Value', 35);
    
    % Chamber Diameter
    uilabel(inputPanel, 'Text', 'Chamber Diameter (mm)', 'Position', [350 30 200 22]);
    dCombEdit = uieditfield(inputPanel, 'numeric', 'Position', [550 30 100 22], 'Value', 130);
    
    % Angle Selection
    uilabel(inputPanel, 'Text', 'Injection Angle (degrees)', 'Position', [10 30 200 22]);
    angleDropdown = uidropdown(inputPanel, 'Items', {'20', '25', '30', '35'}, ...
        'Position', [220 30 100 22], 'Value', '30');
    
    % Max and Min Chamber/Pintle Ratios
    uilabel(inputPanel, 'Text', 'Max Chamber/Pintle Ratio', 'Position', [10 0 200 22]);
    maxRatioEdit = uieditfield(inputPanel, 'numeric', 'Position', [220 0 100 22], 'Value', 5);
    
    uilabel(inputPanel, 'Text', 'Min Chamber/Pintle Ratio', 'Position', [350 0 200 22]);
    minRatioEdit = uieditfield(inputPanel, 'numeric', 'Position', [550 0 100 22], 'Value', 3);
    
    % Compute Button
    computeBtn = uibutton(inputPanel, 'Text', 'Compute', 'Position', [880 20 80 100], ...
        'ButtonPushedFcn', @(btn,event) computeDesign());
    
    % Results Panel
    resultsPanel = uipanel(fig, 'Title', 'Design Results', 'Position', [20 50 960 480]);
    
    % Detailed Results Table
    resultsTable = uitable(resultsPanel, 'Position', [10 10 940 450], ...
        'ColumnName', {'Category', 'Parameter', 'Value', 'Units'});
    
    % Nested function to compute design and display results
    function computeDesign()
        try
            % Collect input values
            dmOx = dmOxEdit.Value;
            OF = ofEdit.Value;
            pInj = pInjEdit.Value;
            pComb = pCombEdit.Value;
            dComb = dCombEdit.Value;
            theta = deg2rad(str2double(angleDropdown.Value));
            maxChamberToPintleRatio = maxRatioEdit.Value;
            minChamberToPintleRatio = minRatioEdit.Value;
            
            % Run full pintle injector design calculation
            [results, dPintle, rPr] = pintleMasterDesign(dmOx, OF, pInj, pComb, dComb, theta, ...
                maxChamberToPintleRatio, minChamberToPintleRatio);
            
            % Prepare comprehensive results data
            data = {
                % Oxidizer Injection Parameters
                'Oxidizer', 'Mass Flow Rate', dmOx, 'kg/s';
                'Oxidizer', 'Injection Pressure', pInj, 'bar'; % 15% already calculated. put something like 50bar
                'Oxidizer', 'Injection Area (aMin)', results.aOx, 'mm²';
                'OxGeometry', 'Pintle Outer Radius (rPost)', results.rPost, 'mm';
                'OxGeometry', 'Open Length (LOpen)', results.LOpen, 'mm';
                'OxGeometry', 'Minimum Length (LMin)', results.LMin, 'mm';
                'OxGeometry', 'Wall Thickness (tWall)', results.twall, 'mm';
                'OxGeometry', 'Pintle Support Radius (rPr)', rPr, 'mm';
                
                % Fuel Injection Parameters
                'Fuel', 'Mass Flow Rate', dmOx/OF, 'kg/s';
                'Fuel', 'Injection Area', results.aFuel, 'mm²';
                'FuelGeometry', 'Annulus Thickness', results.annulusThickness, 'mm';
                
                % Geometric Parameters
                'Geometry', 'Chamber Diameter', dComb, 'mm';
                'Geometry', 'Pintle Diameters', dPintle(2), 'mm';
                'Geometry', 'Max Chamber/Pintle Ratio', maxChamberToPintleRatio, '-';
                'Geometry', 'Min Chamber/Pintle Ratio', minChamberToPintleRatio, '-';
                
             
                
                % Additional Metrics
                'Performance', 'Area Cg (aMin*1.5)', results.aCg, 'mm²';
                'Performance', 'Chamber/Pintle Ratio', results.chamberPintleRatio, '-'
            };
            
            set(resultsTable, 'Data', data);
            
        catch ME
            % Display any errors
            uialert(fig, ME.message, 'Calculation Error');
            disp(ME.stack);
        end
    end

    % Initial computation
    computeDesign();
end

function [results, dPintle, rPr] = pintleMasterDesign(dmOx, OF, pInj, pComb, dComb, theta, maxRatio, minRatio)
    % Oxidizer and Fuel Properties
    oxidizer = N2O;  % Assuming N2O is defined in your local workspace
    
    % Mass flow rates
    dmFuel = dmOx/OF; 
    
    % Injection conditions
    cdOx = 0.65; 
    cdFuelHole = 0.65; 
    
    % Fuel properties
    rhoFuel = 789; 
    
    % Compute pintle diameters based on chamber/pintle ratio
    dPintle = linspace(dComb/maxRatio, dComb/minRatio, 3); %[mm]
    
    % Select first pintle diameter (middle of range)
    rPost = dPintle(2)/2;
    
    % NHNE MODEL calculations (similar to original script)
    hLi = interp1(oxidizer.P, oxidizer.hL, pInj); 
    rhoLi = interp1(oxidizer.P, oxidizer.rhoL, pInj); 
    sL = interp1(oxidizer.P, oxidizer.sL, pInj); 

    % Assuming the outlet pressure (using the choke parameter function)
    fChockParamFun = load("fChockParamFun.mat"); 
    poOx = fChockParamFun.f_chock_param(pInj) .* pInj; 

    % More interpolations and calculations from original script
    hLo = interp1(oxidizer.P, oxidizer.hL, poOx); 
    hVo = interp1(oxidizer.P, oxidizer.hV, poOx); 
    sLo = interp1(oxidizer.P, oxidizer.sL, poOx); 
    sVo = interp1(oxidizer.P, oxidizer.sV, poOx); 

    rhoLo = interp1(oxidizer.P, oxidizer.rhoL, poOx); 
    rhoVo = interp1(oxidizer.P, oxidizer.rhoV, poOx); 

    % Isotropic hypothesis
    so = sL; 

    % Computing the quality
    Xo = (so - sLo)./(sVo - sLo); 

    % Computing the enthalpy
    ho = hLo .* (1 - Xo) + hVo .* Xo; 

    % Computing the outlet density
    rhoo = rhoVo .* rhoLo ./ (rhoLo .* Xo + rhoVo .* (1 - Xo)); 

    % Computing PV vapor pressure
    data.ox.tiOx = 298; % tank temperature
    a = 4.37799; 
    b = 621.077; 
    c = -44.659; 
    PV = 10^(a-(b/data.ox.tiOx+c)); 

    % Computing k parameters
    k = sqrt((pInj - poOx)./(PV - poOx)); 

    % Computing the oxidizer injection area
    aOx = dmOx/( (1-(1/1+k))*cdOx*sqrt(2*rhoLi*(pInj-poOx)*1e5) ...
        + (1/1+k)*cdOx*rhoo*sqrt(2*(hLi-ho)) ); %in m 
    aOx = aOx * 1e6; % Convert aOx to mm²

    % SPI model for fuel
    aFuel = dmFuel/(cdFuelHole*sqrt(2*rhoFuel*(pInj-pComb)*1e5)); 
    aFuel = aFuel*1e6; % [mm²]

    % Design procedure
    rPr_fixed = 5; % minimum pintle root radius for structural integrity 
    ratio_AcgAmin = 1.5; 
    aCg = ratio_AcgAmin * aOx; 

    % Structural calculations
    rPr = rPr_fixed; 
    rCg = sqrt((aCg/pi) + rPr^2); 
    twall = rPost - rCg; 
    
    % Compute channel lengths
    tPost = 0.5; % small chamfer between inner and outer wall of first channel 
    LOpen = ( (rPost-tPost) - ( sqrt((rPost - tPost)^2 - aOx*(sin(theta)/pi) ) ) )/ (sin(theta)*cos(theta));
    LMin = LOpen*cos(theta);
    
    % Second channel radius for fuel
    rSecondChannel = sqrt( (aFuel/pi) + rPost^2); 
    
    % Compute chamber to pintle ratio
    chamberPintleRatio = dComb / (rPost*2);
    
    % Package results
    results = struct(...
        'aOx', aOx, ...
        'aFuel', aFuel, ...
        'aCg', aCg, ...
        'rPost', rPost, ...
        'LOpen', LOpen, ...
        'LMin', LMin, ...
        'twall', twall, ...
        'annulusThickness', rSecondChannel - rPost, ...
        'chamberPintleRatio', chamberPintleRatio ...
    );
end

% Main function to run the GUI
function main()
    pintleInjectorDesignGUI();
end