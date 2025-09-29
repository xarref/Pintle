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

fprintf('Property Table Generator for Pintle Injector Calculator\n');
fprintf('====================================================\n\n');

% Check if CoolProp is available
try
    CoolProp.PropsSI('D', 'T', 298, 'P', 101325, 'Ethanol');
    fprintf('CoolProp detected successfully\n\n');
catch ME
    error('CoolProp not found. Please install CoolProp for MATLAB first.');
end

% Generate ethanol table
fprintf('Generating ethanol property table...\n');
ethanol_table = generate_fluid_table('Ethanol', 250:2.5:350, [1,2,5,10,15,20,25,30,40,50,60,70,80,90,100]);

% Generate N2O table  
fprintf('Generating N2O property table...\n');
n2o_table = generate_fluid_table('N2O', 180:2.5:400, [1,2,5,10,15,20,25,30,40,50,60,70,80,90,100]);

% Save tables
fprintf('\nSaving tables...\n');
save_table(ethanol_table, 'ethanol_properties.json');
save_table(n2o_table, 'n2o_properties.json');

% Summary
fprintf('\n====================================================\n');
fprintf('GENERATION COMPLETE\n');
fprintf('====================================================\n');
fprintf('Ethanol table: %d points, %.1f KB\n', length(ethanol_table.data), get_file_size('ethanol_properties.json')/1024);
fprintf('N2O table:     %d points, %.1f KB\n', length(n2o_table.data), get_file_size('n2o_properties.json')/1024);
fprintf('\nOperating ranges:\n');
fprintf('  Ethanol: %.0f-%.0fK, %.0f-%.0f bar\n', ethanol_table.temperature_range, ethanol_table.pressure_range);
fprintf('  N2O:     %.0f-%.0fK, %.0f-%.0f bar\n', n2o_table.temperature_range, n2o_table.pressure_range);
fprintf('\nFiles ready for web deployment!\n');

end

function table = generate_fluid_table(fluid, temperatures, pressures)
% Generate property table for a single fluid

data = [];
valid_points = 0;
total_points = length(temperatures) * length(pressures);

for T = temperatures
    for P_bar = pressures
        try
            P_pa = P_bar * 1e5; % Convert bar to Pa
            
            % Get density from CoolProp
            density = CoolProp.PropsSI('D', 'T', T, 'P', P_pa, fluid);
            
            % Validate result
            if density > 0 && isfinite(density)
                data(end+1,:) = [T, P_bar, density];
                valid_points = valid_points + 1;
            else
                fprintf('  Invalid density at T=%.1fK, P=%.1fbar: %.3f\n', T, P_bar, density);
            end
            
        catch ME
            fprintf('  Error at T=%.1fK, P=%.1fbar: %s\n', T, P_bar, ME.message);
        end
    end
end

fprintf('  Generated %d/%d valid data points\n', valid_points, total_points);

% Create table structure
table = struct();
table.fluid = lower(fluid);
table.property = 'density';
table.units = 'kg/m3';
table.temperature_range = [min(temperatures), max(temperatures)];
table.temperature_unit = 'K';
table.pressure_range = [min(pressures), max(pressures)];
table.pressure_unit = 'bar';
table.accuracy = '±0.1% (CoolProp reference quality)';
table.data_points = valid_points;
table.data = data; % Format: [T_K, P_bar, rho_kg/m3]

end

function save_table(table, filename)
% Save table to JSON file

% Convert to JSON string
json_str = jsonencode(table);

% Write to file
fid = fopen(filename, 'w');
if fid == -1
    error('Could not open file %s for writing', filename);
end
fprintf(fid, '%s', json_str);
fclose(fid);

% Report file size
file_info = dir(filename);
fprintf('  Saved %s: %.1f KB\n', filename, file_info.bytes/1024);

end

function size_bytes = get_file_size(filename)
% Get file size in bytes

file_info = dir(filename);
if isempty(file_info)
    size_bytes = 0;
else
    size_bytes = file_info.bytes;
end

end
