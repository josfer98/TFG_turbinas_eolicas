% José Antonio Fernández López
% TFG - Generación de energía con una turbina eólica flotante

%% Pruebas de programación

%% Configuración
% Constantes 
    % Densidad del aire
        Ro = 1.225; %Kg/m^3

% Variables
    % Longitud de la pala
        L = 5; %m
    % Número de segmentos en los que se divide la pala
        N = 5;
    % Ángulo inicial de giro
        Theta_1 = 5; %Grados [º]
    % Variación de Theta_i
        Delta_theta = 10; %Grados [º]
    % Velocidad del viento
        u = [1 3 5 8 10 12 15 20]; %Km/h
    % Chord line, de momento sin _i
        c = [2] %m
    % Ángulo de reducción de la Chord line a lo largo de la pala de la
    % turbina eólica
        Phi = 10 %Grados [º]

%% Fórmulas 
    % Lado inicial de la pala
        c_left_pala = c + (L/2) * tan(Phi);
    % Lado final de la pala
        c_right_pala = c + (L/2) * tan(Phi);
    % Área de la pala
        a = ((c_left_pala + c_right_pala) * L) / 2;
    % Fuerza del viento
        F_viento = (1/2) * Ro * a * u^2;
    % Fuerza normal
        F_normal = F_viento * sen(Theta_1);