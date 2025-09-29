function generate_fluid_library()
% GENERATE_FLUID_LIBRARY - Generate comprehensive fluid property library
% 
% This function generates a single comprehensive fluid library JSON file
% containing high-accuracy property data for ethanol and N2O using CoolProp.
%
% Requirements:
%   - CoolProp for MATLAB: https://coolprop.sourceforge.net/coolprop/wrappers/MATLAB/index.html
%
% Output file:
%   - fluid_library.json (~200KB total)

fprintf('Fluid Property Library Generator\n');
fprintf('===============================\n\n');

% Check if CoolProp is available
try
    py.CoolProp.CoolProp.PropsSI('D', 'T', 298, 'P', 101325, 'Ethanol');
    fprintf('CoolProp detected successfully\n\n');
catch ME
    error('CoolProp not found. Please install CoolProp for MATLAB first.');
end

% Define temperature and pressure grids for comprehensive coverage
T_ethanol = 250:2.5:350;  % K, fine resolution for rocket propellant range
T_n2o = 180:2.5:400;      % K, includes supercritical region
P_range = [1, 2, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100]; % bar

% Generate comprehensive fluid library
fprintf('Generating comprehensive fluid library...\n');
fluid_library = struct();

% Library metadata
fluid_library.info = struct();
fluid_library.info.generator = 'MATLAB CoolProp Interface';
fluid_library.info.version = '1.0';
fluid_library.info.date = datestr(now, 'yyyy-mm-dd HH:MM:SS');
fluid_library.info.description = 'High-accuracy fluid property library for rocket propellant calculations';
fluid_library.info.accuracy = 'Reference quality (CoolProp/NIST equivalent)';
fluid_library.info.total_points = 0;

% Generate ethanol data
fprintf('  Processing ethanol...\n');
ethanol_data = generate_fluid_data('Ethanol', T_ethanol, P_range);
fluid_library.ethanol = ethanol_data;
fluid_library.info.total_points = fluid_library.info.total_points + length(ethanol_data.density);

% Generate N2O data  
fprintf('  Processing N2O...\n');
n2o_data = generate_fluid_data('N2O', T_n2o, P_range);
fluid_library.n2o = n2o_data;
fluid_library.info.total_points = fluid_library.info.total_points + length(n2o_data.density);

% Save the complete library
filename = 'fluid_library.json';
fprintf('\nSaving fluid library...\n');
save_json_library(fluid_library, filename);

% Summary
file_info = dir(filename);
fprintf('\n===============================\n');
fprintf('LIBRARY GENERATION COMPLETE\n');
fprintf('===============================\n');
fprintf('Total data points: %d\n', fluid_library.info.total_points);
fprintf('File size: %.1f KB\n', file_info.bytes/1024);
fprintf('Output file: %s\n', filename);
fprintf('\nFluid coverage:\n');
fprintf('  Ethanol: %.0f-%.0fK, %.0f-%.0f bar (%d points)\n', ...
    min(ethanol_data.temperature_grid), max(ethanol_data.temperature_grid), ...
    min(ethanol_data.pressure_grid), max(ethanol_data.pressure_grid), ...
    length(ethanol_data.density));
fprintf('  N2O:     %.0f-%.0fK, %.0f-%.0f bar (%d points)\n', ...
    min(n2o_data.temperature_grid), max(n2o_data.temperature_grid), ...
    min(n2o_data.pressure_grid), max(n2o_data.pressure_grid), ...
    length(n2o_data.density));
fprintf('\nReady for web deployment!\n');

end

function fluid_data = generate_fluid_data(fluid_name, temperatures, pressures)
% Generate comprehensive property data for a single fluid

% Pre-allocate arrays for efficiency
n_points = length(temperatures) * length(pressures);
T_grid = zeros(n_points, 1);
P_grid = zeros(n_points, 1);
rho_data = zeros(n_points, 1);

valid_count = 0;
total_count = 0;

% Generate all T,P combinations
for i = 1:length(temperatures)
    T = temperatures(i);
    for j = 1:length(pressures)
        P_bar = pressures(j);
        P_pa = P_bar * 1e5;
        total_count = total_count + 1;
        
        try
            % Get density from CoolProp
            rho = py.CoolProp.CoolProp.PropsSI('D', 'T', T, 'P', P_pa, fluid_name);
            
            % Validate result
            if rho > 0 && isfinite(rho) && rho < 5000  % Reasonable density limits
                valid_count = valid_count + 1;
                T_grid(valid_count) = T;
                P_grid(valid_count) = P_bar;
                rho_data(valid_count) = rho;
            else
                fprintf('    Invalid density for %s at T=%.1fK, P=%.1fbar: %.3f\n', ...
                    fluid_name, T, P_bar, rho);
            end
            
        catch ME
            fprintf('    Error for %s at T=%.1fK, P=%.1fbar: %s\n', ...
                fluid_name, T, P_bar, ME.message);
        end
    end
end

% Trim arrays to actual valid data
T_grid = T_grid(1:valid_count);
P_grid = P_grid(1:valid_count);
rho_data = rho_data(1:valid_count);

fprintf('    Generated %d/%d valid points for %s\n', valid_count, total_count, fluid_name);

% Create structured fluid data
fluid_data = struct();
fluid_data.name = lower(fluid_name);
fluid_data.temperature_grid = T_grid;
fluid_data.pressure_grid = P_grid;
fluid_data.density = rho_data;
fluid_data.units = struct();
fluid_data.units.temperature = 'K';
fluid_data.units.pressure = 'bar';
fluid_data.units.density = 'kg/m3';
fluid_data.range = struct();
fluid_data.range.temperature = [min(T_grid), max(T_grid)];
fluid_data.range.pressure = [min(P_grid), max(P_grid)];
fluid_data.data_points = valid_count;
fluid_data.accuracy = 'Reference quality (CoolProp)';

end

function save_json_library(library, filename)
% Save the fluid library as a compact JSON file

% Convert to JSON with compact formatting
json_options = struct();
json_options.Compact = true;  % Minimize file size

try
    % Use MATLAB's built-in JSON encoder (R2016b+)
    json_str = jsonencode(library, json_options);
catch
    % Fallback for older MATLAB versions
    fprintf('Warning: Using basic JSON encoding (update MATLAB for better compression)\n');
    json_str = jsonencode(library);
end

% Write to file
fid = fopen(filename, 'w');
if fid == -1
    error('Could not open file %s for writing', filename);
end

fprintf(fid, '%s', json_str);
fclose(fid);

% Report success
file_info = dir(filename);
fprintf('  Saved %s: %.1f KB\n', filename, file_info.bytes/1024);

end
