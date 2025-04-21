function Le = MinThreadLength(pitch, r_uts, h_uts, ring_d, force)
    % Minimum thread length calculation for a ring to prevent thread 
    % from stripping

    % Conversion constants
    in = 25.4;
    lb = 0.22480894;
    psi = 145.037738;

    % Conversions
    n = in / pitch;
    d = ring_d / in;
    force = force * lb;

    % Thread dimensions
    H = sqrt(3)/2 * pitch / in; % mm, height of the thread
    d1 = d - 5 * H / 4; % Minor diameter of external thread
    d2 = d - 3 * H / 4; % Pitch diameter

    % Shear strengths (60% of ultimate tensile strength)
    sh_uts = 0.6 * h_uts * psi; % MPa, shear strength of the housing
    sr_uts = 0.6 * r_uts * psi; % MPa, shear strength of the ring

    % Calculate minimum engagement lengths by solving for Le (pag. 1537
    % Machinist Handbook)
    Le_r = force / (pi * n * d1 * (1 / (2 * n) + 0.57735 * (d2 - d1)) * sr_uts);
    Le_h = force / (pi * n * d * (1 / (2 * n) + 0.57735 * (d - d2)) * sh_uts);

    % The minimum thread length needed to prevent stripping
    Le = max(Le_r, Le_h) * in;
end
