% José Antonio Fernández López
% TFG - Generación de energía con una turbina eólica flotante

%% Pruebas

%% Configuración

    % Densidad del aire
       Ro = 1.225; %Kg/m^3

    % Longitud de la pala
        L = 55; %m

    % Número de segmentos en los que se divide la pala
        N = 5;
        i = 1:N;

    % Ángulo inicial de giro
        theta_1 = 150; %Grados [º]

    % Variación de Theta_i
        Delta_theta = 0.05; %Grados [º]

        %Creamos el ángulo de torsión
        theta_i = zeros(1,N);
    
            for j = 1:N
                if j < 2
                    theta_i(1) = theta_1;
                else
                    theta_i(j) = theta_i(j-1) + Delta_theta;
                end
            end
    
        theta_1 = (theta_1 * pi) / 180; %Rad
        Delta_theta = (Delta_theta * pi) / 180; %Rad
        theta_i = (theta_i .* pi) / 180; %Rad
        
    % Velocidad del viento
        u = 1:1:100; %m/s
        v = 15;
    % Densidad del material de la pala
        CFRP = 14100; %kg/m^3
        GFRP = 15000; %kg/m^3
        GFEpoxi = 17000; %kg/m^3
        dens_pala = [CFRP GFRP GFEpoxi];
    
    % Masa de la pala
        masa_pala = 5000; %Kg

%% Fórmulas para el cálculo inicial y completo de la pala de la turbina eólica.
 
  % Establecemos el buje y la punta
    Buje = 3; %m
    Punta = 1; %m
    L_i = L/N;
  % Se calcula la hipotenusa de borde de fuga
    H_bf = sqrt(((Buje - Punta)^2) + L^2);
  % Ahora el ángulo Phi, con el que decrece la chord line a lo largo de L
    Phi = asin( (Buje - Punta) / H_bf );
    Phi_deg = (Phi * 180) / pi;
  % Variables necesarias para el cálculo de la chord line
    altura_i = (((2*i) -1) * L) / (2*N);
    diagonal_i = (((2*i) -1) * H_bf) / (2*N);
    x_i = sqrt(diagonal_i.^2 - altura_i.^2);
   % Ya se puede obtener la línea de cuerda de cada segmento
    c_i = Buje - x_i;
    
  % Usando algunas fórmulas del desarrollo de Carlos Armenta Deu,
  % referenciado en mi trabajo.
    % Lado inicial de la pala
        c_left_i = c_i + (L_i/2) * tan(Phi);
    % Lado final de la pala
        c_right_i = c_i - (L_i/2) * tan(Phi);
    % Área de la pala
        S_i = ((c_left_i + c_right_i) * L) / 2;

    %Definición del brazo
        cateto_buje = (Buje/2) - (Punta/2);
        R_brazo = sqrt(cateto_buje.^2 + L.^2);
        brazo_i = (((2*i) -1) .* R_brazo) / (2*N);
    % Momento de inercia de un rectángulo
        S_pala = sum(S_i); 
        m_i = (S_i/S_pala) * masa_pala;


        
    % Momento inercia del área de un trapecio
        I_area = (L_i^3).*((c_right_i.^2) + (4.*c_right_i.*c_left_i) + (c_left_i.^2)) ./ (36 .* (c_right_i + c_left_i));
        steiner_theorem = m_i .* (brazo_i.^2);

        I_eje_paralelo = I_area + steiner_theorem;

    % Momento de inercia general
        I = dens_pala(3) .* I_eje_paralelo;



    % Calculamos el tiempo, lo que tarda el viento en atravesar el
    % segmento
        intervalo_tiempo = c_i ./ u(v);
    
%% Cuando solo presenta ángulo de cabeceo
    % Fuerza del viento
        F_viento_i = (1/2) .* Ro .* S_i .* (u(v).^2);
    % Fuerza normal
        F_normal_i = F_viento_i .* sin(theta_1);
    % Momento de torsión
        torque_0 = F_normal_i .* brazo_i;
    % Torque global
        torque_global_0 = sum(torque_0);



    % Aceleración angular
        alpha_ang_0 = torque_global_0 ./ I;
    % Velocidad angular
        Omega_0 = alpha_ang_0 .* intervalo_tiempo;
    % Potencia de la pala
        potencia_0 = torque_global_0 * Omega_0

%% Cuando presenta ángulo de cabeceo y torsión de los segmentos

    % Fuerza normal
        F_normal_i_torsion = F_viento_i .* sin(theta_i);
    % Momento de torsión
        %Se necesita multiplicar por el seno del ángulo de torsión para
        %sacar el torque del área efectiva
        torque_1 = zeros(1,N);
      
            for j = 1:5
                if j < 2
                    torque_1(1) = F_normal_i_torsion(1) .* brazo_i(1);
                else
                       %AQUÍ TMB SE ME PRESENTA UN PROBLEMA, SE ME HACE LA
                       %MEDIDA MUY ENORME, ALGO PASA.
                    torque_1(j) = F_normal_i_torsion(j) .* brazo_i(j) ./ cos(Delta_theta);
                end
            end

    % Torque global
        torque_global_1 = sum(torque_1);

    % Aceleración angular
        alpha_ang_1 = torque_1 ./ I;
    % Velocidad angular
        Omega_1 = alpha_ang_1 .* intervalo_tiempo;
    % Potencia de la pala
        potencia_1 = torque_global_1 * Omega_1


        
 %% Representaciones

    figure('Name','Potencia en cada segmento de la pala con u fija','NumberTitle','off');
        plot(i,potencia_1);
        hold on;
        plot(i,potencia_0);
        xlabel('Intervalo de tiempo (s)');
        ylabel('Potencia (W)');
        legend('Potencia de torsión','Potencia básica')


