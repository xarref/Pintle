%% Main script for the calculation of the minimum thread length 
r_uts = 515; % MPa, ultimate tensile strength of the ring
h_uts = 1200; % MPa, ultimate tensile strength of the housing
ring_d = 48; % mm
pitch = 1.0; % mm
internal_pressure = 6; % MPa
safety_factor = 5;

force = internal_pressure * pi * (ring_d/2)^2 * safety_factor; % N, force on the ring

Le = MinThreadLength(pitch, r_uts, h_uts, ring_d, force);

disp(Le)
