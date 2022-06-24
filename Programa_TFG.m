
%% José Antonio Fernández López
% TFG - Generación de energía con una turbina eólica flotante
%% Variables Globales

    % Número de segmentos en los que se divide la pala
        N = 5;
    % Longitud de la pala
        L = 75; %m
    % Longitud de los segmentos
        L_i = L/N; % m
    % Iteraciones de cálculo iguales al número de segmentos
        i = 1:N;
    % Tamaños del buje y la punta
        BUJE = 3; % m
        PUNTA = 0.6; % m
    % Diametro de la góndola del rotor
        DIAMETRO_GONDOLA = 4; % m
    % Velocidades del viento
        U_VIENTO = 0:1:26;
        % Longitud del vector viento
           M = length(U_VIENTO);
    % Densidad del material de la pala
        CFRP = 1750; %kg/m^3 1.5-2 es el rango bueno
        GFRP = 1050; %kg/m^3 0.91-1.2 es el rango bueno
        GFEPOXY = 1176; %kg/m^3 1159, 1184 y 1186 son los datos, 1176 la mean
        DENS_PALA = [CFRP GFRP GFEPOXY];
    % Densidad del aire
        RHO = 1.225; %Kg/m^3
    % Ancho de la pala en el buje y en la punta
        ANCHO_BUJE = 2; %m
        ANCHO_PUNTA = 0.250; %m

%% Main

    % Llamamos a la función que obtiene el CP mediante regresión polinómica
        CP = coeficiente_potencia;


%  % Contadores para la leyenda
%      contador0 = 1;
%      contador1 = 1;
%
%  figure('Name','Variación ángulo cabeceo 0 - 3.2')
%      for THETA_1_SETUP = 0
% 
%         theta_CL = THETA_1_SETUP;
%         DELTA_THETA = 0.04;
% 
%         [theta_i, THETA_1, DELTA_THETA] = setup_torsion(N, THETA_1_SETUP, DELTA_THETA);
% 
%         [c_left_i, c_right_i, s_i, brazo_i, c_i] = medidas_geometricas(BUJE, PUNTA, L, i, N, L_i);
%             
%         [v_frustum_i, v_frustum_total] = volumen_pala(ANCHO_BUJE, ANCHO_PUNTA, L, i, N, L_i, c_right_i, c_left_i);
%             
%         [I, masa_pala] = momento_inercia(v_frustum_i, s_i, L_i, c_i, c_right_i, c_left_i, DENS_PALA, v_frustum_total, brazo_i);
%         
%         F_viento_i = fuerza_viento(theta_CL, N, M, RHO, L, U_VIENTO);
%          
%         [torque_0, torque_global_0] = torque_cabeceo(F_viento_i, THETA_1, brazo_i);
%         
%         [torque_1, torque_global_1] = torque_torsion(F_viento_i, theta_i, M, N, brazo_i, DELTA_THETA);
% 
%         omega = velocidad_angular(L, DIAMETRO_GONDOLA, BUJE, RHO, U_VIENTO, CP, I);
%         
%         [potencia_0, potencia_1, eta] = potencia_y_eficiencia(omega, torque_global_0, torque_global_1);
% 
%         %eta_1(contador0,:) = eta;
%         
%         
%         [contador0, contador1] = plots(U_VIENTO, potencia_0, potencia_1, contador0, contador1);
%         title('Pot en base a U VIENTO, Variación ángulo cabeceo 0 - 3.2');
% 
%      end
%  % Comandos de la leyenda de la primera figura
%      hold off;
%      legend show;
%      legend('Location','best')

 % Contadores para la leyenda
    contador0 = 0.01;
    contador1 = 0.01;
 figure('Name','Ángulo cabeceo = 2° y ángulo torsión = 0.01°')
    for DELTA_THETA_SETUP = 0.01
        
        THETA_1 = 2;
        theta_CL = THETA_1;
        
        [theta_i, THETA_1, DELTA_THETA] = setup_torsion(N, THETA_1, DELTA_THETA_SETUP);

        [c_left_i, c_right_i, s_i, brazo_i, c_i] = medidas_geometricas(BUJE, PUNTA, L, i, N, L_i);
            
        [v_frustum_i, v_frustum_total] = volumen_pala(ANCHO_BUJE, ANCHO_PUNTA, L, i, N, L_i, c_right_i, c_left_i);
            
        [I, masa_pala] = momento_inercia(v_frustum_i, s_i, L_i, c_i, c_right_i, c_left_i, DENS_PALA, v_frustum_total, brazo_i);
        
        F_viento_i = fuerza_viento(theta_CL, N, M, RHO, L, U_VIENTO);
         
        [torque_0, torque_global_0] = torque_cabeceo(F_viento_i, THETA_1, brazo_i);
        
        [torque_1, torque_global_1] = torque_torsion(F_viento_i, theta_i, M, N, brazo_i, DELTA_THETA);

        omega = velocidad_angular(L, DIAMETRO_GONDOLA, BUJE, RHO, U_VIENTO, CP, I);
        
        [potencia_0, potencia_1, eta] = potencia_y_eficiencia(omega, torque_global_0, torque_global_1, CP);
        
        [contador0, contador1] = plot_torsion(U_VIENTO, potencia_0, potencia_1, contador0, contador1);
        title('Pot en base a U VIENTO, Ángulo cabeceo = 2° y ángulo torsión = 0.01 - 0.04°');

    end
 % Comandos de la leyenda de la segunda figura
    hold off;
    legend show;
    legend('Location','best')
% Para ver si hay algún valor distinto en la eficiencia
%     cuenta=1;
% for prueba=min(U_VIENTO):0.1:max(U_VIENTO)
%     if eta(cuenta) == 1.014
%          disp(cuenta)
%     end
%     cuenta = cuenta + 1;
% end

% Para obtener la eficiencia cuando U_VIENTO tiene determinado valor
% contar = 1;
% for valor_viento=min(U_VIENTO):0.1:max(U_VIENTO)
%          if (valor_viento == 0.8 || valor_viento == 2.4 || valor_viento == 4.3 || valor_viento == 6.7 || ...
%              valor_viento == 9.3 || valor_viento == 12.3 || valor_viento == 15.5 || valor_viento == 18.9)
%             disp(eta(contar))
%          end
%     contar = contar + 1;
% end


% %Prueba para comprobar el CP de mi rotor, se usa uno genérico
% p_rotor = (1/2)*RHO*0.59*pi*(L^2) * 30.^3;
% diametro_rotor = (L * 2) + DIAMETRO_GONDOLA;
% grosor_rotor = BUJE; %m, el grosor depende de la parte mas ancha de la pala
% volumen_rotor = (pi/4) * (diametro_rotor^2) * grosor_rotor;
% p_wind = (1/2)*RHO* volumen_rotor .* 30^3;
% cp = p_rotor ./ p_wind


%  % Contadores para la leyenda
%     contador0 = 1;
%     contador1 = 1;
%  % Cambiamos la longitud de la pala para ver como afecta
%     L = 37.5;
%     L_i = L/N;
% 
% 
%  figure('Name','Hacemos que la pala mida la mitad L=37.5')
%     for THETA_1_SETUP = 0:0.4:3.2
%         theta_CL = THETA_1_SETUP;
%         DELTA_THETA = 0.04;
% 
%         [theta_i, THETA_1, DELTA_THETA] = setup_torsion(N, THETA_1_SETUP, DELTA_THETA);
% 
%         [c_left_i, c_right_i, s_i, brazo_i, c_i] = medidas_geometricas(BUJE, PUNTA, L, i, N, L_i);
%             
%         [v_frustum_i, v_frustum_total] = volumen_pala(ANCHO_BUJE, ANCHO_PUNTA, L, i, N, L_i, c_right_i, c_left_i);
%             
%         [I, masa_pala] = momento_inercia(v_frustum_i, s_i, L_i, c_i, c_right_i, c_left_i, DENS_PALA, v_frustum_total, brazo_i);
%         
%         F_viento_i = fuerza_viento(theta_CL, N, M, RHO, L, U_VIENTO);
%          
%         [torque_0, torque_global_0] = torque_cabeceo(F_viento_i, THETA_1, brazo_i);
%         
%         [torque_1, torque_global_1] = torque_torsion(F_viento_i, theta_i, M, N, brazo_i, DELTA_THETA);
% 
%         omega = velocidad_angular(L, DIAMETRO_GONDOLA, BUJE, RHO, U_VIENTO, CP, I);
%         
%         [potencia_0, potencia_1, eta] = potencia_y_eficiencia(omega, torque_global_0, torque_global_1);
%         
%         [contador0, contador1] = plots(U_VIENTO, potencia_0, potencia_1, contador0, contador1);
%         title('Pot en base a U VIENTO, Hacemos que la pala mida la mitad L=37.5');
%
%         %eta_3(contador0,:) = eta;        
%     end
% 
%  % Comandos de la leyenda de la segunda figura
%     hold off;
%     legend show;
%     legend('Location','best')

 

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
    
    function [c_left_i, c_right_i, s_i, brazo_i, c_i] = medidas_geometricas(BUJE, PUNTA, L, i, N, L_i)
    %% Cálculos de las medidas geométricas de una pala
    
    % Se calcula la hipotenusa de borde de fuga
        h_bf = sqrt(((BUJE - PUNTA)^2) + L^2); % m
    % Ahora el ángulo Phi, con el que decrece la chord line a lo largo de L
        Phi = asin( (BUJE - PUNTA) / h_bf ); % [º]
        % Me servía para conocer cual es el decrecimiento en grados.
        %Phi_deg = (Phi * 180) / pi;
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
    
    function [v_frustum_i, v_frustum_total] = volumen_pala(ANCHO_BUJE, ANCHO_PUNTA, L, i, N, L_i, c_right_i, c_left_i)
    %% Cálculo volumen pala
    
        % Cálculo del volumen del frustum piramidal irregular
            recta_decrecimiento = sqrt((L^2) + (ANCHO_PUNTA - ANCHO_BUJE)^2);
            recta_decrecimiento_i = (i * recta_decrecimiento) / N; % m
        
        % Variables auxiliares para la obtención del área de las bases, z_i
        % para las menores y b_i para las mayores
            z_i = sqrt(recta_decrecimiento_i.^2 - (L_i*i).^2);
            b_i = [0, z_i(1:N-1)];
        
        
        % Ya se puede obtener la línea de cuerda de cada segmento
            ancho_bases_menores = ANCHO_BUJE - z_i; % m
            ancho_bases_mayores = ANCHO_BUJE - b_i; % m
        
        % Con el ancho de las bases y el largo de los segmentos se puede obtener el
        % área de cada uno de los segmentos de la pala
            area_base_menor = ancho_bases_menores .* c_right_i; % m^2
            area_base_mayor = ancho_bases_mayores .* c_left_i; % m^2
        
        % Se calcula el volumen del tronco de pirámide
            v_frustum_i = (L_i/3) .* (area_base_mayor + area_base_menor + sqrt(area_base_menor .* area_base_mayor));
            v_frustum_total = sum(v_frustum_i); % Kg * m^3
    end
    
    function [I, masa_pala] = momento_inercia(v_frustum_i, s_i, L_i, c_i, c_right_i, c_left_i, DENS_PALA, v_frustum_total, brazo_i)
    %%  Momento inercia general de los segmentos de la pala
        
            % Masa de cada segmento de la pala
            s_pala = sum(s_i); % Área de la pala (m)

            % Supuestamente del volumen total de la pala, solo está relleno
            % cerca de un 20%, el resto es aire.
            masa_pala = DENS_PALA(3) * (v_frustum_total*0.2); %Kg
            m_i = (s_i/s_pala) * masa_pala; %Kg de cada segmento
    
            %Ahora se halla la densidad volumétrica
                dens_volumetrica_i = m_i./(v_frustum_i*0.2);
            % Con la dens volumétrica se puede obtener la superficial que
            % se busca
                dens_superficial = dens_volumetrica_i .* c_i;

            % Inercia del área de un trapecio
                I_area = (L_i^3).*((c_right_i.^2) + (4.*c_right_i.*c_left_i) + (c_left_i.^2)) ./ (36 .* (c_right_i + c_left_i));
                I_cm = dens_superficial .* I_area;
               
            % Teorema de Steiner    
                steiner_theorem = m_i .* (brazo_i.^2);
            
            % Momento de inercia general
                I = I_cm + steiner_theorem;
        
    end
        
    function F_viento_i = fuerza_viento(theta_CL, N, M, RHO, L, U_VIENTO)
    %% Fuerza del viento para distintas velocidades

        % Coeficiente de sustentación, mediante réplica de la curva
        % genérica
            if theta_CL >= 180
                theta_CL = theta_CL - 180;
            end

            if theta_CL <= 12 
                C_L = ((1.1/12) * (theta_CL)) + 0.3; % Reynolds
            elseif theta_CL <= 17 
                C_L = ((1.18-1.4)/(17-12) * (theta_CL - 12)) + 1.4; % Reynolds
            elseif theta_CL <= 45
                C_L = ((1.79-1.18)/(45-17) * (theta_CL - 17)) + 1.18; % Reynolds
            elseif theta_CL <= 90
                C_L = ((-1.79)/(90-45) * (theta_CL - 45)) + 1.79; % Reynolds
            elseif theta_CL <= 135
                C_L = ((-0.9)/(135-90) * (theta_CL - 90)); % Reynolds
            elseif theta_CL <= 180-17
                C_L = ((-1.18+1.79)/(163-135) * (theta_CL - 135)) - 1.79; % Reynolds
            elseif theta_CL <= 180-12
                C_L = ((-1.79+1.18)/(168-163) * (theta_CL - 163)) - 1.18; % Reynolds
            elseif theta_CL <= 180
                C_L = ((0+1.4)/(180-168) * (theta_CL - 168)) - 1.4; % Reynolds            
            end


        F_viento_i = zeros(N,1);
        for j = 1:M
            F_viento_i(j) = (1/2) .* RHO .* ((pi/4)*((2*L).^2)) .* (U_VIENTO(j).^2) * C_L;
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

    function omega = velocidad_angular(L, DIAMETRO_GONDOLA, BUJE, RHO, U_VIENTO, CP, I)
%% velocidad_angular


    % Se va a calular la energía cinética generada en la pala, para ello se
    % define el disco que definen las palas al realizar una vuelta completa en el que se mide la masa de aire que llega al rotor
    % Diámetro del rotor completo
        diametro_rotor = (L * 2) + DIAMETRO_GONDOLA; % m
    
    % Grosor del disco
        grosor_rotor = BUJE; %m, el grosor depende de la parte mas ancha de la pala
    
    % Volumen del disco
        volumen_rotor = (pi/4) * (diametro_rotor^2) * grosor_rotor; %m^3
    
    % Masa de aire que atraviesa el volumen del disco
        M_aire = volumen_rotor * RHO; %Kg

    % Energía cinética
        %e_cinetica = 1/2 * M_aire * U_VIENTO^2;
    % Involucramos el coeficiente de potencia con la energía cinética
         %e_cinetica_aprovechada = e_cinetica * CP;
    % Ahora necesitamos la energía cinética de rotación
         %e_cinetica_rotacion = 1/2 * I * omega^2;

    % Se compara y despeja con la energía cinética aprovechada para obtener
    % la velocidad angular omega
        omega = sqrt(M_aire./I) .* U_VIENTO.';
        omega = sqrt(CP.') .* omega;
end
    
    function [potencia_0, potencia_1, eta] = potencia_y_eficiencia(omega, torque_global_0, torque_global_1, CP)
    %% Desarrollo potencia y eficiencia
        % Potencia de la pala
            potencia_0 = (torque_global_0 .* CP.') .* omega;
            potencia_1 = (torque_global_1 .* CP.') .* omega;

            %potencia_0 = (torque_global_0) .* omega;
            %potencia_1 = (torque_global_1) .* omega;
            
            potencia_0 = sum(potencia_0,2);
            potencia_1 = sum(potencia_1,2);
            
        % Se calcula el % de mejora o empeoramiento mediante la torsión de los
        % segmentos de la pala
            eta = potencia_1.' ./ potencia_0.';
    end
    
    function  [contador0, contador1] = plot_torsion(U_VIENTO, potencia_0, potencia_1, contador0, contador1)
%% Representaciones

%Potencia obtenida dependiendo de la velocidad del viento
    x = U_VIENTO;
    y0 = potencia_0.'/1e6;
    y1 = potencia_1.'/1e6;
    % Necesito hacer un contador externo para la leyenda ya que el show de la
    % leyenda está fuera del propio bucle de potencias que se están mostrando.
    % Si no se hace esto, MATLAB ignora la extra entries de la leyenda y solo
    % muestra la primera.

    %leyenda0 = ['Potencia SIN torsión ', num2str(contador0)]; contador0 = 1 + contador0;
    leyenda1 = ['Potencia CON torsión ', num2str(contador1),'°']; contador1 = 0.01 + contador1;
    leyenda0 = ['Potencia SIN torsión '];
    
    % Para que solo haga plot 1 vez de la potencia sin torsión
    if contador1 == 0.02
    plot(x, y0, 'DisplayName', leyenda0,'LineStyle','-'); hold on;
    end

    
    hold on;plot(x, y1, 'DisplayName', leyenda1,'LineStyle','--');
    
    xlabel('Velocidad del viento (m/s)');
    ylabel('Potencia (MW)');
    xlim([min(U_VIENTO) max(U_VIENTO)])

    end

    function  [contador0, contador1] = plots(U_VIENTO, potencia_0, potencia_1, contador0, contador1)
%% Representaciones

%Potencia obtenida dependiendo de la velocidad del viento
    x = U_VIENTO;
    y0 = potencia_0.'/1e6;
    y1 = potencia_1.'/1e6;
    % Necesito hacer un contador externo para la leyenda ya que el show de la
    % leyenda está fuera del propio bucle de potencias que se están mostrando.
    % Si no se hace esto, MATLAB ignora la extra entries de la leyenda y solo
    % muestra la primera.

    leyenda0 = ['Potencia SIN torsión ', num2str(contador0)]; contador0 = 1 + contador0;
    leyenda1 = ['Potencia CON torsión ', num2str(contador1)]; contador1 = 1 + contador1;
    
   
    plot(x, y0, 'DisplayName', leyenda0);
    hold on;plot(x, y1, 'DisplayName', leyenda1);
    
    xlabel('Velocidad del viento (m/s)');
    ylabel('Potencia (MW)');
    xlim([min(U_VIENTO) max(U_VIENTO)])

    end

    function CP = coeficiente_potencia
%% Coeficiente de potencia
% Se necesita el coeficient de potencia en base a la función del
% viento, lo obtenemos de uno genérico.
    y = [0.01	0.01	0.01	0.02	0.05	0.28	0.4	0.48	0.48	0.45	0.4	0.35	0.3	0.23	0.19	0.15	0.12	0.1	0.09	0.08	0.07	0.06  0.05	0.04	0.035	0.03	0.027];
    x = 0:1:26;
    
    % Se realiza una regresión polinomial para obtener la curva del CP en
    % base a uno genérico usado en los NACA
    p = polyfit(x,y,14);
    
    CP = zeros (1,26);
    var_aux = 0;
    for vel_viento = 1:1:27
        if var_aux < 3 || var_aux == 26
            CP(1,vel_viento) = 0;
        else
            CP(1,vel_viento) = p(1)  .* (vel_viento.^14) + p(2)  .* (vel_viento.^13) + p(3)  .* (vel_viento.^12)...
                + p(4)  .* (vel_viento.^11) + p(5)  .* (vel_viento.^10) + p(6)  .* (vel_viento.^9)...
                + p(7)  .* (vel_viento.^8)  + p(8)  .* (vel_viento.^7)  + p(9)  .* (vel_viento.^6)...
                + p(10) .* (vel_viento.^5)  + p(11) .* (vel_viento.^4)  + p(12) .* (vel_viento.^3)...
                + p(13) .* (vel_viento.^2)  + p(14) .* (vel_viento)     + p(15);
        end
        var_aux = var_aux + 1;
    end
end