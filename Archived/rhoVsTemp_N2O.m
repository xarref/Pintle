% MATLAB script to plot density vs temperature for N2O using CoolProp

% Temperature range in Celsius
temperature_celsius = -15:0.1:35; % Range from -15 to 30 with a step of 0.1

%fluid
fluid = 'N2O';

% Convert temperature to Kelvin (K = C + 273.15)
temperature_kelvin = temperature_celsius + 273.15;

% Initialize arrays to store density values
density_quality_0 = zeros(size(temperature_kelvin));
density_quality_1 = zeros(size(temperature_kelvin));

% Define pressure range (from 1 bar to 51 bar in increments of 5 bar)
pressure_range = 1e5:5e5:70e5; % Pressure in Pa

density_pressure = zeros(length(pressure_range), length(temperature_kelvin));

% Use CoolProp to calculate density for N2O (Nitrous Oxide)
for i = 1:length(temperature_kelvin)
    % Density calculated using quality (Q=0 and Q=1)
    density_quality_0(i) = py.CoolProp.CoolProp.PropsSI('D', 'T', temperature_kelvin(i), 'Q', 0, fluid); % Saturated liquid density
    density_quality_1(i) = py.CoolProp.CoolProp.PropsSI('D', 'T', temperature_kelvin(i), 'Q', 1, fluid); % Saturated vapor density
    
    for j = 1:length(pressure_range)
        % Density calculated using pressure
        density_pressure(j, i) = py.CoolProp.CoolProp.PropsSI('D', 'T', temperature_kelvin(i), 'P', pressure_range(j), fluid);
    end
end

% Plot density vs temperature for all pressures
figure;
plot(temperature_celsius, density_quality_0, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Density (Q=0)');
hold on;
plot(temperature_celsius, density_quality_1, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Density (Q=1)');

colors = lines(length(pressure_range)); % Generate distinct colors for each pressure
for j = 1:length(pressure_range)
    plot(temperature_celsius, density_pressure(j, :), '--', 'LineWidth', 1.5, 'Color', colors(j, :), 'DisplayName', sprintf('P = %d bar', pressure_range(j)/1e5));
end
hold off;





grid on;
title('Density vs Temperature for N2O');
xlabel('Temperature (°C)');
ylabel('Density (kg/m^3)');

% Customize plot
legend('show', 'Location', 'Best');
