function rocket_simulation()
    % Initialize parameters
    params = initialize_rocket_parameters();
    
    T = readtable('thrust_curve.csv', ...
        'VariableNamingRule','preserve', ...
        'TreatAsEmpty',{'NA','','–'});  % treat blanks or dashes as empty

% Identify which columns are numeric
% Assume the first two numeric columns are time and thrust
varTypes = varfun(@class, T, 'OutputFormat','cell');
isNumCol = strcmp(varTypes,'double') | strcmp(varTypes,'single');

numCols = find(isNumCol);
if numel(numCols) < 2
    error('Thrust CSV must contain at least two numeric columns (time, thrust).');
end

timeCol   = numCols(1);
thrustCol = numCols(2);

% Extract data and remove any rows with NaN/Inf
rawTime   = T{:, timeCol};
rawThrust = T{:, thrustCol};
maskValid = isfinite(rawTime) & isfinite(rawThrust);

thrust_time = rawTime(maskValid);
thrust_vals = rawThrust(maskValid);

% Sort by time (interp1 requires monotonic X)
[thrust_time, sortIdx] = sort(thrust_time);
thrust_vals = thrust_vals(sortIdx);

% Store in params
params.thrust_time = thrust_time;
params.thrust_vals = thrust_vals;
    % User-configurable launch angle (degrees)
    params.launch_angle = 90;  
    
    % Set up simulation
    tspan = [0 params.max_time];
    vx0 = params.initial_speed * cosd(params.launch_angle);
    vy0 = params.initial_speed * sind(params.launch_angle);
    initial_conditions = [0; vx0; 0; vy0];  % [x, vx, y, vy]
    
    % Configure ODE solver options
    options = odeset(...
        'Events', @(t,y) event_functions(t,y,params), ...
        'RelTol',1e-6,'AbsTol',1e-8);
    
    % Run simulation
    [t, y, te, ye, ie] = ode45(@(t,y) rocket_dynamics(t,y,params), ...
                               tspan, initial_conditions, options);
    
        % ------------ New Code to Sample and Save State Data ------------

    % Desired sampling rate (Hz) and interval
    fs = 10;            % samples per second
    dt = 1/fs;          % seconds between samples
    t_sample = 0:dt:t(end);

    % Interpolate each state onto the uniform grid
    x_sample   = interp1(t, y(:,1), t_sample);
    vx_sample  = interp1(t, y(:,2), t_sample);
    alt_sample = interp1(t, y(:,3), t_sample);
    vy_sample  = interp1(t, y(:,4), t_sample);
    
    % Compute flight path angle (degrees) at each sample
    angle_sample = atan2d(vy_sample, vx_sample);
    
    % Package into a MATLAB timeseries
    ts_x   = timeseries(x_sample,   t_sample, 'Name','x');
    ts_vx  = timeseries(vx_sample,  t_sample, 'Name','vx');
    ts_alt = timeseries(alt_sample, t_sample, 'Name','altitude');
    ts_vy  = timeseries(vy_sample,  t_sample, 'Name','vy');
    ts_angle = timeseries(angle_sample,   t_sample, 'Name','flight_angle');

    % Combine into a single timedataset or save to .mat for Simulink import
    save('rocket_states.mat','ts_x','ts_vx','ts_alt','ts_vy','ts_angle');

    % Or write to CSV for external import
    T = table(t_sample(:), x_sample(:), vx_sample(:), alt_sample(:), vy_sample(:), angle_sample(:), ...
          'VariableNames', {'time','x','vx','altitude','vy','flight_angle'});

    writetable(T,'rocket_states.csv');

    % ---------------- End of Sampling & Saving --------------------
    
    % Real-time visualization
    visualize_flight(t, y, te, ye, ie, params);

end

function dydt = rocket_dynamics(t, y, params)
    x = y(1); vx = y(2); alt = y(3); vy = y(4);
    
    % Determine flight phase
    phase = determine_flight_phase(t, alt, vy, params);
    
    % Compute thrust, mass, drag coeff
    [thrust, mass, cd, ref_area] = calculate_forces(t, phase, params);
    
    % Atmosphere & gravity
    [rho, g] = atmosphere_model(alt);
    
    % Wind disturbance (small sinusoidal gust)
    wind = params.wind_amp * sin(2*pi*params.wind_freq*t);
    % in rocket_dynamics, after calculating phase:
    if phase >= 4  % main chute deployed
        effective_wind = params.wind_amp_parachute * sin(2*pi*params.wind_freq*t); % much smaller
    else
        effective_wind = params.wind_amp * sin(2*pi*params.wind_freq*t);
    end

    % Relative velocity
    vrel_x = vx - effective_wind;
    vrel_y = vy;
    vrel = sqrt(vrel_x^2 + vrel_y^2);
    
    % Drag

    Fd = 0.5 * rho * vrel^2 * ref_area * cd;

    if vrel>0
        drag_x = Fd * (vrel_x/vrel);
        drag_y = Fd * (vrel_y/vrel);
    else
        drag_x = 0; drag_y = 0;
    end
    
    % Equations of motion
    dxdt  = vx;
    dvxdt = (thrust*cosd(params.launch_angle) - drag_x) / mass;
    dydt  = vy;
    dvydt = (thrust*sind(params.launch_angle) - mass*g - drag_y) / mass;
    
    dydt = [dxdt; dvxdt; dydt; dvydt];
end

function phase = determine_flight_phase(t, alt, vy, params)
    persistent reached_apogee deployed_drogue deployed_main
    if isempty(reached_apogee),    reached_apogee=false; end
    if isempty(deployed_drogue),    deployed_drogue=false; end
    if isempty(deployed_main),      deployed_main=false; end
    
    if t <= params.burn_time
        phase = 1;                              % powered ascent
    elseif ~reached_apogee && vy<=0
        phase = 2; reached_apogee=true;         % apogee
    elseif reached_apogee && ~deployed_drogue
        phase = 3; deployed_drogue=true;        % drogue deployed at apogee
    elseif reached_apogee && deployed_drogue && ~deployed_main && alt<=params.main_deploy_alt
        phase = 4; deployed_main=true;          % main deployed on descent
    elseif deployed_main
        phase = 5;                              % main-parachute descent
    else
        phase = 2;                              % coasting
    end
end

function [thrust, mass, cd, ref_area] = calculate_forces(t, phase, params)
    switch phase
        case 1  % powered ascent
            thrust = interp1(params.thrust_time, params.thrust_vals, t, 'linear', 0);
            mass   = params.init_mass - params.mass_flow*t;
            cd     = params.cd_before_main;
            ref_area = params.ref_area_before_main;
        case 2  % coasting to apogee
            thrust = 0;
            mass   = params.dry_mass;
            cd     = params.cd_before_main;
            ref_area = params.ref_area_before_main;
        case 3  % drogue descent
            thrust = 0;
            mass   = params.dry_mass;
            cd     = params.cd_before_main;
            ref_area = params.ref_area_before_main;

        case {4,5}  % main descent (main chute deployed)
            thrust = 0;
            mass   = params.dry_mass;
            cd     = params.cd_after_main;
            ref_area = params.ref_area_after_main;
    end
end


function [value,isterm,direction] = event_functions(t,y,params)
    alt = y(3); vy = y(4);
    % 1 = apogee (vy=0 downwards)
    value(1)=vy;      isterm(1)=0; direction(1)=-1;
    % 2 = drogue at apogee (trigger immediately after apogee)
    value(2)=vy + eps;isterm(2)=0; direction(2)=-1;
    % 3 = main parachute altitude
    value(3)=alt-params.main_deploy_alt; isterm(3)=0; direction(3)=-1;
    % 4 = ground impact
    value(4)=alt;    isterm(4)=1; direction(4)=-1;
end

function visualize_flight(t, y, te, ye, ie, params)
    figure('Position', [100 100 1200 800]);
    % Subplots setup
    ax1 = subplot(2,3,[1,4]); hold(ax1,'on'); grid(ax1,'on');
    xlabel(ax1, 'Downrange (m)'); ylabel(ax1, 'Altitude (m)');
    title(ax1, 'Trajectory');

    ax2 = subplot(2,3,2); grid(ax2,'on');
    xlabel(ax2, 'Time (s)'); ylabel(ax2, 'Altitude (m)');
    title(ax2, 'Altitude vs Time');

    ax3 = subplot(2,3,3); grid(ax3,'on');
    speed = sqrt(y(:,2).^2 + y(:,4).^2);
    plot(ax3, t, speed, 'k-'); xlabel(ax3, 'Time (s)');
    ylabel(ax3, 'Speed (m/s)'); title(ax3, 'Speed vs Time');

    % Remove/comment out PHASE plot
    % ax4 = subplot(2,3,5); grid(ax4,'on');
    % phases = zeros(size(t));
    % for i = 1:numel(t)
    %     phases(i) = determine_flight_phase(t(i), y(i,3), y(i,4), params);
    % end
    % stairs(ax4, t, phases, 'm-', 'LineWidth', 1.5);
    % xlabel(ax4, 'Time (s)'); ylabel(ax4, 'Phase');
    % title(ax4, 'Flight Phase'); ylim(ax4, [0.5 5.5]); yticks(ax4, 1:5);

    % Add Wind disturbance plot
    ax5 = subplot(2,3,5); grid(ax5,'on');
    wind = params.wind_amp * sin(2*pi*params.wind_freq*t);
    plot(ax5, t, wind, 'b-'); xlabel(ax5, 'Time (s)');
    ylabel(ax5, 'Wind (m/s)'); title(ax5, 'Wind Disturbance');

    ax6 = subplot(2,3,6); grid(ax6, 'on');
    plot(ax6, t, y(:,4), 'g-');
    xlabel(ax6, 'Time (s)');
    ylabel(ax6, 'Vertical Velocity (m/s)');
    title(ax6, 'Vertical Velocity vs Time');

    % Real-time animation
    for i = 1:length(t)
        % Trajectory plot
        subplot(2,3,[1,4]);
        plot(y(1:i,1), y(1:i,3),'b-','LineWidth',2);
        h = plot(y(i,1), y(i,3),'ro','MarkerFaceColor','r');
        % --- Altitude vs Time plot inside the loop ---
        subplot(2,3,2);
        cla; % Clear axis for fresh plot
        plot(t(1:i), y(1:i,3), 'r-', 'LineWidth', 1.2);
        xlabel('Time (s)');
        ylabel('Altitude (m)');
        title('Altitude vs Time');
        grid on;
        pause(0.01);
        delete(h);
        drawnow;
    end
end



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


function [rho,g] = atmosphere_model(alt)
    if alt<=11000
        T=288.15-0.0065*alt;
        P=101325*(T/288.15)^5.256;
        rho=P/(287*T);
    else
        rho=0.3639*exp(-0.000157*alt);
    end
    R=6.371e6; g=9.81*(R/(R+alt))^2;
end
