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
        %Theta_1 = 5; %Grados [º]
        Theta_1 = (5 * pi) / 180; %Rad
    % Variación de Theta_i
        %Delta_theta = 10; %Grados [º]
        Delta_theta = (10 * pi) / 180; %Rad
    % Velocidad del viento
        u = [5 8 10 12]; %m/s
    % Chord line, de momento sin _i
        c = [2]; %m
    % Ángulo de reducción de la Chord line a lo largo de la pala de la
    % turbina eólica
        %Phi = 20; Grados [º]
        Phi = (20 * pi) / 180; %Rad


    % Tiempo, se puede usar como algo concreto o como array para
    % representaciones
        tiempo = 0; % ¿Qué debemos asignarle?

%% Fórmulas para el cálculo inicial y completo de la pala de la turbina eólica.
  % Usando algunas fórmulas del desarrollo de Carlos Armenta Deu,
  % referenciado en mi trabajo.
    % Lado inicial de la pala
        c_left_pala = c + (L/2) * tan(Phi);
    % Lado final de la pala
        c_right_pala = c - (L/2) * tan(Phi);
    % Área de la pala
        a = ((c_left_pala + c_right_pala) * L) / 2;
    
    % Fuerza del viento
        F_viento = (1/2) * Ro * a * u.^2;
    % Fuerza normal
        F_normal = F_viento * sin(Theta_1);
    % Brazo := distancia entre el inicio de la pala y el centro de giro.
        %En este caso voy a suponer un valor hasta conocer como se debe
        %calcular, el valor debe ser aprox ~0.5-0.6 la pala
        brazo_inicial = 0.55 * L;
    % Momento de torsión
        torque_0 = F_normal * brazo_inicial
    
    % Momento de inercia
        % I = masa_pala * (brazo_inicial^2)
    % Aceleración angular
        % alpha_ang = torque_0 / I;
    % Velocidad angular
        % Omega = alpha_ang * tiempo;
    
    % Potencia de la pala
        % potencia_0 = torque_0 * Omega;


%% Representaciones
    

    %Desconozco si esta es una correcta representación, en uno de los
    %estudios se representaba con xlabel de tip speed.
    
    figure('Name','Potencia vs. Velocidad del viento','NumberTitle','off');
        plot(u,potencia_0);
        xlabel('Velocidad del viento (m/s)');
        ylabel('Potencia de la pala (W)');



