function params = initialize_rocket_parameters()
    params.init_mass           = 25.0;    % kg
    params.dry_mass            = 20.0;    % kg
    params.mass_flow           = 1.25;    % kg/s
    params.ref_area            = 0.05;    % m^2
    params.main_chute_area = 1.0; % m^2 (or something much larger than 0.05)
    params.burn_time           = 4.0;     % s
    params.max_thrust          = 500;     % N
    params.initial_speed       = 0;       % m/s
    params.cd_body             = 0.15;
    params.cd_drogue           = 1;
    params.cd_main             = 2.5;
    params.main_deploy_alt     = 100;     % m
    params.max_time            = 200;     % s
    params.wind_amp            = 3;       % m/s
    params.wind_freq           = 0.1;     % Hz
    params.launch_angle        = 90;      % deg (default)
    params.cd_before_main = 0.15;          % Cd – before main deployment (body phases)
    params.cd_after_main  = 2.5;           % Cd – after main deployment (parachute)
    params.ref_area_before_main = 0.05;    % m^2 – body area
    params.ref_area_after_main  = 0.5;     % m^2 – approx. parachute canopy area
    params.wind_amp_parachute = 0.01; % much gentler drift (or set to average wind)

end