%% PINTLE INJECTOR DESIGN SCRIPT
% FUNCTION NAME: 
% pintleMasterV4
%
% DESCRIPTION:
%   Script for analyzing and designing pintle injector configurations 
%   with two different geometrical arrangements (Oxidizer Internal/Fuel External 
%   and Fuel Central/Oxidizer External)
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

clc; clear;

%% SET
% Constant parameters decided a priori

% Injection and Combustion Parameters

SET.pInj = 60 * 0.95;  % Injection pressure, assuming 5% pressure drop in lines [bar]
SET.pComb = 35;        % Combustion chamber pressure [bar]
SET.dComb = 72;        % Combustion chamber diameter [mm] (rigenerativo) 
SET.OF = 2.9;          % Oxidizer to Fuel ratio [-]

% Fluid Properties
SET.oxidizer = 'N2O';   % Oxidizer fluid name
SET.fuel = 'Ethanol';   % Fuel fluid name
SET.dmOx = 1.15;        % Mass flow rate Oxidizer [kg/s]
SET.cdOx = 0.65;        % Discharge coefficient for oxidizer [-]
SET.tankTemperature = 298; % Tank temperature [K]

% Geometric Constraints
SET.maxChamberToPintleRatio = 5; 
SET.minChamberToPintleRatio = 3; 

% Fuel Specific Parameters
SET.dmFuel = SET.dmOx/SET.OF;  % Mass flow rate Fuel [kg/s]
SET.cdFuelHole = 0.611;  % Discharge coefficient for fuel holes [-]

% Design Parameters
SET.rPr_fixed = 4;     % Minimum pintle root radius for structural integrity [mm]
SET.ratio_AcgAmin = 1.25;  % Ratio of chamber to minimum injection area [-]
SET.tPost = 1.5;       % Small chamfer between inner and outer wall of first channel [mm]



function result = PropSI(output, input1, value1, input2, value2, fluid)
    result = py.CoolProp.CoolProp.PropsSI(output, input1, value1, input2, value2, fluid);
end

% Arrays to store results for both configurations
results = cell(2,1);

% Run both configurations
for configuration = 1:2
    % Modify pintle diameter range based on configuration allowing for
    % different geometries
    if configuration == 1 % Ox Internal, Fuel External
        dPintle = linspace(SET.dComb/SET.maxChamberToPintleRatio, SET.dComb/SET.minChamberToPintleRatio, 5); %[mm]
    else % Fuel Central, Ox External
        dPintle = linspace(SET.dComb/SET.minChamberToPintleRatio, SET.dComb/SET.maxChamberToPintleRatio, 5); %[mm]
    end

    % Fuel density computation
    rhoFuel = PropSI('D', 'P', (SET.pInj-6)*1e5, 'T', 490, SET.fuel); %assuming 6 bar pressure drop
  
    %% NHNE MODEL
    % Interpolating liquid and vapor properties prior injection
    % Using PropSI for interpolation
    hLi = PropSI('H', 'P', SET.pInj*1e5, 'Q', 0, SET.oxidizer); % Liquid enthalpy at injection pressure
    rhoLi = PropSI('D', 'P', SET.pInj*1e5, 'Q', 0, SET.oxidizer); % Liquid density at injection pressure
    sL = PropSI('S', 'P', SET.pInj*1e5, 'Q', 0, SET.oxidizer); % Liquid entropy at injection pressure

    %%% Computing the outlet conditions
    % Assuming the outlet pressure
    fChockParamFun = load("fChockParamFun.mat"); 
    poOx = fChockParamFun.f_chock_param(SET.pInj) .* SET.pInj; 

    % Interpolating liquid and vapor properties after injection
    hLo = PropSI('H', 'P', poOx*1e5, 'Q', 0, SET.oxidizer); 
    hVo = PropSI('H', 'P', poOx*1e5, 'Q', 1, SET.oxidizer); 
    sLo = PropSI('S', 'P', poOx*1e5, 'Q', 0, SET.oxidizer); 
    sVo = PropSI('S', 'P', poOx*1e5, 'Q', 1, SET.oxidizer); 

    rhoLo = PropSI('D', 'P', poOx*1e5, 'Q', 0, SET.oxidizer); 
    rhoVo = PropSI('D', 'P', poOx*1e5, 'Q', 1, SET.oxidizer); 

    % Isotropic hypothesis
    so = sL; 
    % Computing the quality
    Xo = (so - sLo)./(sVo - sLo); 
    % Computing the enthalpy
    ho = hLo .* (1 - Xo) + hVo .* Xo; 
    % Computing the outlet density
    rhoo = rhoVo .* rhoLo ./ (rhoLo .* Xo + rhoVo .* (1 - Xo)); 
    % Computing PV vapor pressure
    PV = PropSI('P', 'T', SET.tankTemperature, 'Q', 1, SET.oxidizer) / 1e5; % Convert Pa to bar
    % Computing k parameters
    k = sqrt((SET.pInj - poOx)./(PV - poOx)); 

    %% Calculate Areas (NHNE and SPI)
    % Oxidizer injection area (internal) (NHNE)
    aOx = SET.dmOx/( (1-(1/1+k))*SET.cdOx*sqrt(2*rhoLi*(SET.pInj-poOx)*1e5) ...
        + (1/1+k)*SET.cdOx*rhoo*sqrt(2*(hLi-ho)) ); %in m 
    aOx = aOx * 1e6; % Convert aOx to mm 

    % Fuel injection area (external) (SPI)
    aFuel = SET.dmFuel/(SET.cdFuelHole*sqrt(2*rhoFuel*(SET.pInj-SET.pComb)*1e5)); 
    aFuel = aFuel*1e6; % [mm] 

    %%% Computing the injection areas (configuration-dependent)
    if configuration == 1 % Ox Internal, Fuel External
        aMin = aOx;  % Set aMin to aOx for first configuration where ox is central
    else % Fuel Central, Ox External
        aMin = aFuel; % Set aMin to aFuel for second configuration
    end

    %% DESIGN PROCEDURE
    % Ratio of Acg to Amin (design parameter)
    aCg = SET.ratio_AcgAmin*aMin; 

    % Initialize arrays to store results
    rPost = dPintle/2; % Convert diameter to radius 
    rCg = zeros(size(rPost)); 
    rPr = zeros(size(rPost)); 
    twall_first_channel = zeros(size(rPost)); 

    % Design angle options
    theta = [deg2rad(10), deg2rad(20), deg2rad(30), deg2rad(35)]; 

    % Calculate Rcg, Rpr, and wall thickness for each rPost
    for i = 1:length(rPost)
        % Fix Rpr to the structural requirement
        rPr(i) = SET.rPr_fixed; 

        % Calculate Rcg using the inverse formula
        rCg(i) = sqrt((aCg/pi) + rPr(i)^2); 

        % Calculate wall thickness
        twall_first_channel(i) = rPost(i) - rCg(i); 

        % Validate results
        if twall_first_channel(i) < 1
            warning(['Wall thickness for rPost = ', num2str(rPost(i)), ' mm is dangerously small.']); 
        end
    end

    % Compute channel open length
    LOpen = zeros(length(rPost), length(theta)); 
    LMin = zeros(length(rPost), length(theta)); 

    for i = 1:length(rPost)
        for j = 1:length(theta)
            % Compute channel open length based on new geometry
            LOpen(i,j) = ( (rPost(i)-SET.tPost) - ( sqrt((rPost(i) - SET.tPost)^2 - aMin*(sin(theta(j))/pi) ) ) )/ (sin(theta(j))*cos(theta(j))); 
            LMin(i,j) = LOpen(i,j)*cos(theta(j)); 
        end
    end

    % Compute second channel radius (configuration-dependent)
    if configuration == 1 % Ox Internal, Fuel External
        rSecondChannel = zeros(length(rPost),1); 
        for i = 1:length(rPost)
            rSecondChannel(i) = sqrt( (aFuel/pi) + rPost(i)^2); 
        end
    else % Fuel Central, Ox External
        rSecondChannel = zeros(length(rPost),1); 
        for i = 1:length(rPost)
            rSecondChannel(i) = sqrt( (aOx/pi) + rPost(i)^2); 
        end
    end

    % Store results for this configuration
    results{configuration}.configuration = configuration;
    results{configuration}.aMin = aMin;
    results{configuration}.aFuel = (configuration == 1) * aFuel + (configuration ~= 1) * aOx;
    results{configuration}.aCg = aCg;
    results{configuration}.rPost = rPost;
    results{configuration}.rPr = rPr;
    results{configuration}.rSecondChannel = rSecondChannel;
    results{configuration}.twall_first_channel = twall_first_channel;
    results{configuration}.LOpen = LOpen;
    results{configuration}.LMin = LMin;
    results{configuration}.theta = theta;
end

%% Comparative Output
fprintf('Comparative Analysis of Pintle Injector Configurations\n');
fprintf('===================================================\n\n');

for configuration = 1:2
    if configuration == 1
        configName = 'Oxidizer Internal, Fuel External';
    else
        configName = 'Fuel Central, Oxidizer External';
    end
    
    fprintf('Configuration %d: %s\n', configuration, configName);
    fprintf('---------------------------------------------\n');
    
    % Display aMin, aFuel, aCg
    fprintf('Computed aMin (primary injection area): %.4f mm²\n', results{configuration}.aMin); 
    fprintf('Computed secondary injection area: %.4f mm²\n', results{configuration}.aFuel); 
    fprintf('Computed aCg: %.3f mm² (Amin x %.4f) \n\n', results{configuration}.aCg, SET.ratio_AcgAmin); 

    % Display rPr for all rPost values
    fprintf('Computed radius of the internal pintle for the outer radius:\n');
    fprintf('------------------------------------------------\n');
    fprintf('| rPost (mm) |    rPr (mm)  | chamb/pintle (-) |\n');
    fprintf('------------------------------------------------\n');
    for i = 1:length(results{configuration}.rPost)
        fprintf('| %8.2f   | %12.4f | %12.4f     |\n', ...
            results{configuration}.rPost(i), ...
            results{configuration}.rPr(i), ...
            SET.dComb./(results{configuration}.rPost(i)*2) );
    end
    fprintf('------------------------------------------------\n\n');

    % Display other parameters for each combination of rPost and theta
    fprintf('Parameters for each combination of rPost and theta:\n');
    fprintf('-------------------------------------------------------------------------------------------\n');
    fprintf('| rPost (mm) | Theta (deg) | LOpen (mm) | LMin (mm)  | Twall (mm) | Annulus thickness (mm) |\n'); 
    fprintf('-------------------------------------------------------------------------------------------\n');
    for i = 1:length(results{configuration}.rPost)
        for j = 1:length(results{configuration}.theta)
            fprintf('| %8.2f   | %8.2f    | %10.4f | %10.4f | %10.4f | %16.4f       |\n', ...
                results{configuration}.rPost(i), ...
                rad2deg(results{configuration}.theta(j)), ...
                results{configuration}.LOpen(i, j), ...
                results{configuration}.LMin(i,j), ...
                results{configuration}.twall_first_channel(i), ...
                results{configuration}.rSecondChannel(i)-results{configuration}.rPost(i)); 
        end
    end
    fprintf('-------------------------------------------------------------------------------------------\n\n');
end

% %% Grafical Rapresentation
% 
% % Ask user for configuration choice: 1 for Oxidizer Internal, 2 for Fuel Internal
% scelta = input('Configuration (1 for Oxidizer Internal, 2 for Fuel Internal): ');
% 
% if scelta == 1 || scelta == 2
%     % Loop over each pintle radius (rPost) for plotting
%     for i = 1:length(rPost)
%         figure;
%         sgtitle(['Pintle Injector Geometry for rPost = ', num2str(rPost(i)), ' mm']);
%         hold on;
%         axis equal;
%         grid on;
% 
%         % Loop over each design angle (theta) for plotting
%         for j = 1:length(theta)
%             % Current pintle radius and angle
%             rPost_current = rPost(i); 
%             theta_current = theta(j);
% 
%             % Dimensions for the current configuration
%             pintle_radius = rPr(i); % Using fixed pintle root radius from the setup
%             first_channel_inner_width = rCg(i); % Inner width for first channel
%             first_channel_outer_width = rPost_current; % Outer width for first channel
%             pintle_tip_length = (first_channel_outer_width - pintle_radius) * tan(theta_current); % Tip length calculation
%             second_channel_radius = rSecondChannel(i); % Radius of the second channel
%             injector_length = 50; % Injector length assumption
%             t_thickness_tip = 2.5; % Tip thickness
% 
%             % Color assignment based on configuration
%             if scelta == 1
%                 color1 = 'b'; % Oxidizer color (blue)
%                 color2 = 'r'; % Fuel color (red)
%             else
%                 color1 = 'r'; % Fuel color (red)
%                 color2 = 'b'; % Oxidizer color (blue)
%             end
% 
%             % Outer and inner bounds of the first channel
%             first_inner_top = first_channel_inner_width;
%             first_inner_bottom = -first_channel_inner_width;
%             first_outer_top = first_channel_outer_width;
%             first_outer_bottom = -first_channel_outer_width;
%             tWall_top = first_inner_top + SET.tPost; % Wall thickness at the top of the first channel
%             tWall_bottom = first_inner_bottom - SET.tPost; % Wall thickness at the bottom of the first channel
% 
%             % Bounds for the second channel
%             second_channel_top = second_channel_radius;
%             second_channel_bottom = -second_channel_radius;
% 
%             % Splash base for the angled tip
%             splash_base_width = 2 * first_channel_outer_width;
%             splash_top = splash_base_width / 2;
%             splash_bottom = -splash_base_width / 2;
% 
%             % Plotting the geometry for the current configuration and angle
%             subplot(2, 3, j); % Subplot for each angle
%             hold on;
%             axis equal;
%             grid on;
% 
%             % Plot the pintle body (cylindrical part)
%             line([-injector_length, 0], [pintle_radius, pintle_radius], 'Color', 'k', 'LineWidth', 2);
%             line([-injector_length, 0], [-pintle_radius, -pintle_radius], 'Color', 'k', 'LineWidth', 2);
%             line([0, 0], [-pintle_radius, pintle_radius], 'Color', 'k', 'LineWidth', 2);
% 
%             % Plot the angled pintle tip
%             line([0, pintle_tip_length], [pintle_radius, first_outer_top], 'Color', 'k', 'LineWidth', 2);
%             line([0, pintle_tip_length], [-pintle_radius, first_outer_bottom], 'Color', 'k', 'LineWidth', 2);
%             line([pintle_tip_length, t_thickness_tip], [splash_top, splash_top], 'Color', 'k', 'LineWidth', 2);
%             line([pintle_tip_length, t_thickness_tip], [-splash_top, -splash_top], 'Color', 'k', 'LineWidth', 2);
%             line([t_thickness_tip, t_thickness_tip], [splash_bottom, splash_top], 'Color', 'k', 'LineWidth', 2);
% 
%             % Plot the first channel (inner and outer bounds)
%             line([-injector_length, -5], [first_inner_top, first_inner_top], 'Color', color1, 'LineWidth', 2);
%             line([-injector_length, -5], [first_inner_bottom, first_inner_bottom], 'Color', color1, 'LineWidth', 2);
%             line([-injector_length, pintle_tip_length - LOpen(i, j)], [first_outer_top, first_outer_top], 'Color', color1, 'LineWidth', 2);
%             line([-injector_length, pintle_tip_length - LOpen(i, j)], [first_outer_bottom, first_outer_bottom], 'Color', color1, 'LineWidth', 2);
%             line([-5, pintle_tip_length - LOpen(i, j)], [first_inner_top, tWall_top], 'Color', color1, 'LineWidth', 2);
%             line([pintle_tip_length - LOpen(i, j), pintle_tip_length - LOpen(i, j)], [tWall_top, first_outer_top], 'Color', color1, 'LineWidth', 2);
%             line([-5, pintle_tip_length - LOpen(i, j)], [first_inner_bottom, tWall_bottom], 'Color', color1, 'LineWidth', 2);
%             line([pintle_tip_length - LOpen(i, j), pintle_tip_length - LOpen(i, j)], [tWall_bottom, first_outer_bottom], 'Color', color1, 'LineWidth', 2);
% 
%             % Plot the second channel
%             line([-injector_length, pintle_tip_length - LOpen(i, j)], [second_channel_top, second_channel_top], 'Color', color2, 'LineWidth', 2);
%             line([-injector_length, pintle_tip_length - LOpen(i, j)], [second_channel_bottom, second_channel_bottom], 'Color', color2, 'LineWidth', 2);
% 
%             % Add title with specific angle
%             title(['\theta = ', num2str(rad2deg(theta_current)), '^\circ']);
% 
%             % Adjust plot labels
%             xlabel('Length (mm)');
%             ylabel('Height (mm)');
% 
%             hold off;
%         end
%     end
% else
%     % Display an error message if an invalid configuration is selected
%     disp('Invalid choice. Use 1 for Oxidizer Internal or 2 for Fuel Internal.');
% end



%% Third configuration - Fuel internal with holes
% Fuel Hole Configuration Analysis
fprintf('-----------------------------------------------------------\n');
fprintf('Configuration 3: Fuel Internal with Holes and N2O Annulus \n\n\n');
fprintf('Computed fuel injection area: %.4f mm²\n', aFuel); 

% Hole diameter options (meters)
dFuelHoleOptions = [0.3, 0.4, 0.5, 0.6, 0.8] * 1e-3;
rowOptions = [1, 2, 3, 4, 5];

% Prepare storage for results
resultsMatrix = cell(length(dFuelHoleOptions), length(rowOptions), length(rPost));

% Iterate through hole diameters
for i = 1:length(dFuelHoleOptions)
    % Current hole diameter
    dHole = dFuelHoleOptions(i);

    % Compute hole area
    holeArea = pi * (dHole/2)^2;


    % Compute exact number of holes required and then round to the closest
    % integer
    exactHoles = SET.dmFuel / (SET.cdFuelHole * holeArea * sqrt(2*rhoFuel*(SET.pInj - SET.pComb)*1e5));
    maxHoles = round(exactHoles);


    % Analyze distribution across different row configurations
    for j = 1:length(rowOptions)
        nRows = rowOptions(j);

        % Compute holes per row (rounded up)
        holesPerRow = round(maxHoles / nRows);

        % Calculate total holes
        totalHoles = nRows * holesPerRow;

        % Calculate actual mass flow rate
        actualMassFlow = totalHoles * holeArea * SET.cdFuelHole * sqrt(2*rhoFuel*(SET.pInj - SET.pComb)*1e5);

        % Calculate angular spacing between holes
        angularSpacing = 360 / holesPerRow;  % Degrees

        % Convert angular spacing to radians (to calculate arcDistance)
        angularSpacingRad = deg2rad(angularSpacing) - 2*dHole;

        % Iterate through rPost values for each row
        for k = 1:length(rPost)
            % Calculate arc distance for the current rPost value
            arcDistance = rPost(k) * angularSpacingRad;

            % Store results, including the current rPost value
            resultsMatrix{i,j,k} = struct(...
                'holeDiameter', dHole * 1000, ... % mm
                'holeArea', holeArea * 1e6, ... % mm²
                'nRows', nRows, ...
                'holesPerRow', holesPerRow, ...
                'angularSpacing', angularSpacing, ...
                'totalHoles', totalHoles, ...
                'actualMassFlow', actualMassFlow, ...
                'percentError', abs(actualMassFlow - SET.dmFuel) / SET.dmFuel * 100, ...
                'arcDistance', arcDistance, ...
                'rPost', rPost(k) ... % Add rPost value to the structure
            );
        end
    end
end

% Display results in a table format
fprintf('Fuel Discrete Hole Configuration Analysis\n\n');
fprintf('---------------------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf('| Hole Dia. | Hole Area  | Rows | Holes/Row   | Total Holes  | Mass Flow (kg/s) | Error (%%) || Angular Delta | Arc Distance (mm)  | rPost (mm)    |\n');
fprintf('|-------------------------------------------------------------------------------------------------------------------------------------------------|\n');
for i = 1:length(dFuelHoleOptions)
    for j = 1:length(rowOptions)
        for k = 1:length(rPost)
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
                res.rPost ...
            );
        end
    end
    fprintf('|-------------------------------------------------------------------------------------------------------------------------------------------------|\n');
end


%% Performance



    % Calculate velocities using Bernoulli's equation
    v_radial = sqrt(2 * (SET.pInj - SET.pComb) / rhoFuel);
    v_annular = sqrt(2 * (SET.pInj - SET.pComb) / rhoo);
    angle = 45;

    % Calculate TMR
    TMR = ( SET.dmFuel* v_radial) / (SET.dmOx * v_annular);

