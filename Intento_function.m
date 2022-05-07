% José Antonio Fernández López
% TFG - Generación de energía con una turbina eólica flotante
%% Variables Globales

    % Número de segmentos en los que se divide la pala
        N = 5;
    % Longitud de la pala
        L = 34; %m
    % Longitud de los segmentos
        L_i = L/N; % m
    % Iteraciones de cálculo iguales al número de segmentos
        i = 1:N;
    % Ángulo inicial de giro
        %THETA_1 = 1; %Grados [º]
    % Variación de Theta_i para la torsión
        %DELTA_THETA = 0.03; %Grados [º]
    % Tamañós del buje y la punta
        BUJE = 3; % m
        PUNTA = 0.6; % m
    % Velocidades del viento
        U_VIENTO = 1:0.5:20;
        % Longitud del vector viento
           M = length(U_VIENTO);
    % Tiempo de análisis del sistema
        TIEMPO_ANALISIS = 60; %segundos
    % Densidad del material de la pala
        CFRP = 1410; %kg/m^3
        GFRP = 1500; %kg/m^3
        GFEPOXI = 1700; %kg/m^3
        DENS_PALA = [CFRP GFRP GFEPOXI];
    % Densidad del aire
        RHO = 1.225; %Kg/m^3

%% Main

     figure(1)
         for THETA_1_SETUP = 1:0.2:2

            DELTA_THETA = 0.03;

            [theta_i, THETA_1, DELTA_THETA] = setup_torsion(N, THETA_1_SETUP, DELTA_THETA);

            [c_left_i, c_right_i, s_i, brazo_i] = medidas_geometricas(BUJE, PUNTA, L, i, N, L_i);
                
            [v_frustum_i, v_frustum_total] = volumen_pala(L, i, N, L_i, c_right_i, c_left_i);
                
            [I, masa_pala] = momento_inercia(v_frustum_i, s_i, L_i, L, c_right_i, c_left_i, DENS_PALA, v_frustum_total, brazo_i);
            
            F_viento_i = fuerza_viento(N, M, RHO, s_i, U_VIENTO);
             
            [torque_0, torque_global_0] = torque_cabeceo(F_viento_i, THETA_1, brazo_i);
            
            [torque_1, torque_global_1] = torque_torsion(F_viento_i, theta_i, M, N, brazo_i, DELTA_THETA);
            
            [potencia_0, potencia_1] = potencia_y_eficiencia(torque_0, I, torque_1, TIEMPO_ANALISIS, torque_global_0, torque_global_1);
            
            plots(U_VIENTO, potencia_0, potencia_1);
            
            clear -regexp \d
         end
    

    figure(2)
        for DELTA_THETA_SETUP = 0.01:0.02:0.06

            THETA_1 = 1;

            [theta_i, THETA_1, DELTA_THETA] = setup_torsion(N, THETA_1, DELTA_THETA_SETUP);

            [c_left_i, c_right_i, s_i, brazo_i] = medidas_geometricas(BUJE, PUNTA, L, i, N, L_i);
                
            [v_frustum_i, v_frustum_total] = volumen_pala(L, i, N, L_i, c_right_i, c_left_i);
                
            [I, masa_pala] = momento_inercia(v_frustum_i, s_i, L_i, L, c_right_i, c_left_i, DENS_PALA, v_frustum_total, brazo_i);
            
            F_viento_i = fuerza_viento(N, M, RHO, s_i, U_VIENTO);
             
            [torque_0, torque_global_0] = torque_cabeceo(F_viento_i, THETA_1, brazo_i);
            
            [torque_1, torque_global_1] = torque_torsion(F_viento_i, theta_i, M, N, brazo_i, DELTA_THETA);
            
            [potencia_0, potencia_1] = potencia_y_eficiencia(torque_0, I, torque_1, TIEMPO_ANALISIS, torque_global_0, torque_global_1);
            
            plots(U_VIENTO, potencia_0, potencia_1);
            
            clear -regexp \d
        end
    


%% Funciones


    function [theta_i, THETA_1, DELTA_THETA] = setup_torsion(N, THETA_1, DELTA_THETA)
    %% Setup ángulo torsión
    
    %Creación del ángulo de torsión
    theta_i = zeros(1,N);
    for j = 1:N
        if j < 2
            theta_i(1) = THETA_1;
        else
            theta_i(j) = theta_i(j-1) + DELTA_THETA;
        end
    end
    
    THETA_1 =     (THETA_1 * pi)     / 180; %Rad
    DELTA_THETA = (DELTA_THETA * pi) / 180; %Rad
    theta_i =     (theta_i .* pi)    / 180; %Rad
    
    end
    
    function [c_left_i, c_right_i, s_i, brazo_i] = medidas_geometricas(BUJE, PUNTA, L, i, N, L_i)
    %% Cálculos de las medidas geométricas de una pala
    
    % Se calcula la hipotenusa de borde de fuga
    h_bf = sqrt(((BUJE - PUNTA)^2) + L^2); % m
    % Ahora el ángulo Phi, con el que decrece la chord line a lo largo de L
    Phi = asin( (BUJE - PUNTA) / h_bf ); % [º]
    Phi_deg = (Phi * 180) / pi;
    % Variables necesarias para el cálculo de la chord line
    altura_i = (((2*i)-1) * L) / (2*N);
    diagonal_i = (((2*i) -1) * h_bf) / (2*N);
    x_i = sqrt(diagonal_i.^2 - altura_i.^2);
    % Ya se puede obtener la línea de cuerda de cada segmento
    c_i = BUJE - x_i; % m
    
    % Usando algunas fórmulas del desarrollo de Carlos Armenta Deu,
    % referenciado en mi trabajo.
    
    % Lado inicial de la pala
    c_left_i = c_i + (L_i/2) * tan(Phi); % m 
    % Lado final de la pala
    c_right_i = c_i - (L_i/2) * tan(Phi); % m
    % Área de cada segmento de la pala
    s_i = ((c_left_i + c_right_i) / 2) * L_i; % m^2
    
    % Definición del brazo
    cateto_buje = (BUJE/2) - (PUNTA/2);
    r_brazo = sqrt(cateto_buje.^2 + L.^2);
    brazo_i = (((2*i) -1) .* r_brazo) / (2*N); %m
    end
    
    function [v_frustum_i, v_frustum_total] = volumen_pala(L, i, N, L_i, c_right_i, c_left_i)
    %% Cálculo volumen pala
    
    % Cálculo del volumen del frustum piramidal irregular
    ancho_buje = 2; %m
    ancho_punta = 0.250; %m
    
    recta_decrecimiento = sqrt((L^2) + (ancho_punta - ancho_buje)^2);
    recta_decrecimiento_i = (i * recta_decrecimiento) / N; % m
    
    % Variables auxiliares para la obtención del área de las bases, z_i
    % para las menores y b_i para las mayores
    z_i = sqrt(recta_decrecimiento_i.^2 - (L_i*i).^2);
    b_i = [0, z_i(1:N-1)];
    
    
    % Ya se puede obtener la línea de cuerda de cada segmento
    ancho_bases_menores = ancho_buje - z_i; % m
    ancho_bases_mayores = ancho_buje - b_i; % m
    
    % Con el ancho de las bases y el largo de los segmentos se puede obtener el
    % área de cada uno de los segmentos de la pala
    area_base_menor = ancho_bases_menores .* c_right_i; % m^2
    area_base_mayor = ancho_bases_mayores .* c_left_i; % m^2
    
    % Se calcula el volumen del tronco de pirámide mediante las fórmulas del
    % papiro de Moscú
    v_frustum_i = (L_i/3) .* (area_base_mayor + area_base_menor + sqrt(area_base_menor .* area_base_mayor));
    v_frustum_total = sum(v_frustum_i); % Kg * m^3
    % V_frustum = (L/3) * (6 + 0.15 + sqrt(6*0.15)) 
    %Esto es del volumen completo, pero da 80.45 y lo de arriba 80.73, tengo que revisarlo.
    end
    
    function [I, masa_pala] = momento_inercia(v_frustum_i, s_i, L_i, L, c_right_i, c_left_i, DENS_PALA, v_frustum_total, brazo_i)
    %%  Momento inercia general de los segmentos de la pala
    
    % Una vez se obtiene el volumen de la figura, se puede calcular el espesor,
    % https://journals.pan.pl/Content/109465/PDF/AME_125441.pdf
    espesor_i = (v_frustum_i) ./ s_i; %m
    espesor_medio = sum( (espesor_i .* L_i) ./ L ); %m
    
    % Inercia del área de un trapecio
    I_area = (L_i^3).*((c_right_i.^2) + (4.*c_right_i.*c_left_i) + (c_left_i.^2)) ./ (36 .* (c_right_i + c_left_i));
    I_general = espesor_medio .* I_area;
    
    
    % Masa de cada segmento de la pala
    S_pala = sum(s_i); % Área de la pala (m)
    masa_pala = DENS_PALA(1) * (v_frustum_total*0.2); %Kg
    m_i = (s_i/S_pala) * masa_pala; %Kg de cada segmento
    
    % Teorema de Steiner    
    steiner_theorem = m_i .* (brazo_i.^2);
    
    % Momento de inercia general
    I = I_general + steiner_theorem;
    
    
    % Se calcula el tiempo, lo que tarda el viento para diferentes velocidades en atravesar el
    % segmento
    %       intervalo_tiempo = c_i ./ u(v);
    %         intervalo_tiempo = zeros(M,N);
    %         for j = 1:M
    %             for j2 = 1:N
    %                 intervalo_tiempo(j,j2) = c_i(j2) ./ U_VIENTO(j);
    %             end
    %         end
    end
    
    function F_viento_i = fuerza_viento(N, M , RHO, s_i, U_VIENTO)
    %% Fuerza del viento para distintas velocidades
    %F_viento_i = (1/2) .* RHO .* s_i .* (u(v).^2);
    F_viento_i = zeros(M,N);
    for j = 1:M
        for j2 = 1:N
            F_viento_i(j,j2) = (1/2) .* RHO .* s_i(j2) .* U_VIENTO(j);
        end 
    end
    end
    
    function [torque_0, torque_global_0] = torque_cabeceo(F_viento_i, THETA_1, brazo_i)
    %% Torque para ángulo de cabeceo
    % Fuerza normal
    F_normal_i = F_viento_i .* sin(THETA_1);
    % Momento de torsión
    torque_0 = F_normal_i .* brazo_i;
    % Torque global
    torque_global_0 = sum(torque_0,2);
    end
    
    function [torque_1, torque_global_1] = torque_torsion(F_viento_i, theta_i, M, N, brazo_i, DELTA_THETA_prueba)
    %% Torque para ángulo de cabeceo y torsión
    
    % Fuerza normal
    F_normal_i_torsion = F_viento_i .* sin(theta_i);
    
    % Momento de torsión
    torque_1 = zeros(M,N);
    for j = 1:M
        for j2 = 1:N
            if j2 < 2
                torque_1(j,j2) = F_normal_i_torsion(j,j2) .* brazo_i(j2);
            else
                torque_1(j,j2) = F_normal_i_torsion(j,j2) .* brazo_i(j2) .* cos(DELTA_THETA_prueba);
            end
        end
    end
    
    % Torque global
    torque_global_1 = sum(torque_1,2);
    end
    
    function [potencia_0, potencia_1] = potencia_y_eficiencia(torque_0, I, torque_1, TIEMPO_ANALISIS, torque_global_0, torque_global_1)
    %% Desarrollo potencia y eficiencia
    % Aceleración angular
    alpha_ang_0 = torque_0 ./ I;
    alpha_ang_1 = torque_1 ./ I;
    % Velocidad angular
    Omega_0 = TIEMPO_ANALISIS .* alpha_ang_0;
    Omega_1 = TIEMPO_ANALISIS .* alpha_ang_1;
    % Recordar calcular la otra forma!!!!!!!!!!!!!!!!!!!!!!!
    % Potencia de la pala
    potencia_0 = torque_global_0 .* Omega_0;
    potencia_1 = torque_global_1 .* Omega_1;
    
    potencia_0 = sum(potencia_0,2);
    potencia_1 = sum(potencia_1,2);
    
    % Se calcula el % de mejora o empeoramiento mediante la torsión de los
    % segmentos de la pala
    eta = potencia_1 ./ potencia_0;
    end
    
    function plots(U_VIENTO, potencia_0, potencia_1)
%% Representaciones

%Potencia obtenida dependiendo de la velocidad del viento
x = U_VIENTO;
y0 = potencia_0.';
y1 = potencia_1.';

plot(x,y0);hold on;plot(x,y1);

title('Potencia obtenida en 60 segundos en base a la velocidad del viento');
xlabel('Velocidad del viento (m/s)');
ylabel('Potencia (W)');
legend('Potencia SIN torsión 1','Potencia CON torsión 1', ...
    'Potencia SIN torsión 2','Potencia CON torsión 2', ...
    'Potencia SIN torsión 3','Potencia CON torsión 3', ...
    'Potencia SIN torsión 4','Potencia CON torsión 4', ...
    'Potencia SIN torsión 5','Potencia CON torsión 5', ...
    'Potencia SIN torsión 6','Potencia CON torsión 6');
end