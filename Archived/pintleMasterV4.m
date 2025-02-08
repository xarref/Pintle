
%% v4


clc; clear
function result = PropSI(output, input1, value1, input2, value2, fluid)
    result = py.CoolProp.CoolProp.PropsSI(output, input1, value1, input2, value2, fluid);
end

%%% GENERAL DATA
pInj = 60 * 0.85; % assuming 15% pressure drop in the lines 
pComb = 35; 
dComb = 100; % rigenerativo 
OF = 3.8; 

% Ox
oxidizer = 'N2O'; % Change to string for PropSI
fuel = 'Ethanol'; %Change to string for PropSI
dmOx = 1.5; % originariamente 1.3 
cdOx = 0.65; 
tankTemperature = 298; % tank 
maxChamberToPintleRatio = 5; 
minChamberToPintleRatio = 3; 

% Arrays to store results for both configurations
results = cell(2,1);

% Run both configurations
for configuration = 1:2
    % Modify pintle diameter range based on configuration allowing for
    % different geometries
    if configuration == 1 % Ox Internal, Fuel External
        dPintle = linspace(dComb/maxChamberToPintleRatio, dComb/minChamberToPintleRatio, 5); %[mm]
    else % Fuel Central, Ox External
        dPintle = linspace(dComb/minChamberToPintleRatio, dComb/maxChamberToPintleRatio, 5); %[mm]
    end

    % Fuel
    dmFuel = dmOx/OF; % originariamente 0.25 
    cdFuelHole = 0.6; 
    %rhoFuel = 789; %T_amb P_amb
    rhoFuelPropSI = PropSI('D', 'P', pInj*1e5*0.8, 'T', 450, fuel); %density of ethanol at 450K, assuming 20% pressure drop in channels (temperature 
    % is within the safe margin of the specific heat of ethanol)

    %% NHNE MODEL

    % Interpolating liquid and vapor properties prior injection
    % Using PropSI for interpolation
    hLi = PropSI('H', 'P', pInj*1e5, 'Q', 0, oxidizer); % Liquid enthalpy at injection pressure
    rhoLi = PropSI('D', 'P', pInj*1e5, 'Q', 0, oxidizer); % Liquid density at injection pressure
    sL = PropSI('S', 'P', pInj*1e5, 'Q', 0, oxidizer); % Liquid entropy at injection pressure

    %%% Computing the outlet conditions

    % Assuming the outlet pressure
    fChockParamFun = load("fChockParamFun.mat"); 
    poOx = fChockParamFun.f_chock_param(pInj) .* pInj; 

    % Interpolating liquid and vapor properties after injection
    hLo = PropSI('H', 'P', poOx*1e5, 'Q', 0, oxidizer); 
    hVo = PropSI('H', 'P', poOx*1e5, 'Q', 1, oxidizer); 
    sLo = PropSI('S', 'P', poOx*1e5, 'Q', 0, oxidizer); 
    sVo = PropSI('S', 'P', poOx*1e5, 'Q', 1, oxidizer); 

    rhoLo = PropSI('D', 'P', poOx*1e5, 'Q', 0, oxidizer); 
    rhoVo = PropSI('D', 'P', poOx*1e5, 'Q', 1, oxidizer); 

    % Isotropic hypothesis
    so = sL; 
    % Computing the qualityao
    Xo = (so - sLo)./(sVo - sLo); 
    % Computing the enthalpy
    ho = hLo .* (1 - Xo) + hVo .* Xo; 
    % Computing the outlet density
    rhoo = rhoVo .* rhoLo ./ (rhoLo .* Xo + rhoVo .* (1 - Xo)); 
    % Computing PV vapor pressure -> Replace with PropSI vapor pressure
    PV = PropSI('P', 'T', tankTemperature, 'Q', 1, oxidizer) / 1e5; % Convert Pa to bar
    % Computing k parameters
    k = sqrt((pInj - poOx)./(PV - poOx)); 

    %% Calculate Areas (NHNE and SPI)
        % Oxidizer injection area (internal) (NHNE)
        aOx = dmOx/( (1-(1/1+k))*cdOx*sqrt(2*rhoLi*(pInj-poOx)*1e5) ...
        + (1/1+k)*cdOx*rhoo*sqrt(2*(hLi-ho)) ); %in m 
        aOx = aOx * 1e6; % Convert aOx to mm 

        % Fuel injection area (external) (SPI)
        aFuel = dmFuel/(cdFuelHole*sqrt(2*rhoFuelPropSI*(pInj-pComb)*1e5)); 
        aFuel = aFuel*1e6; % [mm] 

    %%% Computing the injection areas (configuration-dependent)
    if configuration == 1 % Ox Internal, Fuel External
        aMin = aOx;  % Set aMin to aOx for first configuration where ox is central
    else % Fuel Central, Ox External
        aMin = aFuel; % Set aMin to aFuel for second configuration
    end

    %% DESIGN PROCEDURE

    % Structural requirement for pintle rigidity
    rPr_fixed = 4; % 5mm minimum pintle root radius for structural integrity 

    % Ratio of Acg to Amin (design parameter)
    ratio_AcgAmin = 1.5; 

    % Recalculate Amin based on the original area calculations
    aCg = ratio_AcgAmin*aMin; 

    % Initialize arrays to store results
    rPost = dPintle/2; % Convert diameter to radius 
    rCg = zeros(size(rPost)); 
    rPr = zeros(size(rPost)); 
    twall_first_channel = zeros(size(rPost)); 

    % Design angle options
    theta = [deg2rad(10), deg2rad(20), deg2rad(30), deg2rad(35)]; 

    % Design parameters
    tPost = 1.5; % small chamfer between inner and outer wall of first channel 

    % Calculate Rcg, Rpr, and wall thickness for each rPost
    for i = 1:length(rPost)
        % Fix Rpr to the structural requirement
        rPr(i) = rPr_fixed; 

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
            LOpen(i,j) = ( (rPost(i)-tPost) - ( sqrt((rPost(i) - tPost)^2 - aMin*(sin(theta(j))/pi) ) ) )/ (sin(theta(j))*cos(theta(j))); 
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
                    %When configuration is 1: (1 * aFuel) + (0 * aOx) = aFuel
                    %When configuration is 2: (0 * aFuel) + (1 * aOx) = aOx
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
    fprintf('Computed aCg: %.4f mm²\n\n', results{configuration}.aCg); 

    % Display rPr for all rPost values
    fprintf('Computed radius of the internal pintle for the outer radius:\n');
    fprintf('------------------------------------------------\n');
    fprintf('| rPost (mm) |    rPr (mm)  | chamb/pintle (-) |\n');
    fprintf('------------------------------------------------\n');
    for i = 1:length(results{configuration}.rPost)
        fprintf('| %8.2f   | %12.4f | %12.4f     |\n', ...
            results{configuration}.rPost(i), ...
            results{configuration}.rPr(i), ...
            dComb./(results{configuration}.rPost(i)*2) );
    end
    fprintf('------------------------------------------------\n\n');

    % Display other parameters for each combination of rPost and theta
    fprintf('Parameters for each combination of rPost and theta:\n');
    fprintf('-------------------------------------------------------------------------------------------\n');
    fprintf('| rPost (mm) | Theta (deg) | LOpen (mm) | LMin (mm)  | Twall (mm) | Annulus thickness (mm) |\n'); 
    fprintf('-------------------------------------------------------------------------------------------\n');
    for i = 1:length(results{configuration}.rPost)
        fprintf('|            |             |            |            |            |                        |\n')
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