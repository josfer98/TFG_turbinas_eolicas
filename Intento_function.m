% José Antonio Fernández López
% TFG - Generación de energía con una turbina eólica flotante
%% Configuración de la función

    % Longitud de la pala
        L = 34; %m
    % Número de segmentos en los que se divide la pala
        N = 5;
    % Ángulo inicial de giro
        %theta_1 = 1; %Grados [º]
    % Variación de Theta_i para la torsión
        %Delta_theta = 0.03; %Grados [º]
    % Masa de la pala
        masa_pala = 5200; %Kg
    % Tamañós del buje y la punta
        Buje = 3; %metros
        Punta = 0.6; %metros
    % Velocidades del viento
        u_viento = 1:0.5:20;
    % Tiempo de análisis del sistema
        tiempo_analisis = 60; %segundos

figure(1)
Delta_theta = 0.03;
for theta_1 = 1:0.2:2
calculo_potencias(N, L, theta_1, Delta_theta, Buje, Punta, masa_pala, u_viento, tiempo_analisis);
end

figure(2)
theta_1 = 1;
for Delta_theta = 0.01:0.01:0.06
calculo_potencias(N, L, theta_1, Delta_theta, Buje, Punta, masa_pala, u_viento, tiempo_analisis);
end




function calculo_potencias(N, L, theta_1, Delta_theta, Buje, Punta, masa_pala, u_viento, tiempo_analisis)
%% Setup

% Densidad del aire
Ro = 1.225; %Kg/m^3
% Iteraciones de cálculo iguales al número de segmentos
i = 1:N;
% Longitud de los segmentos
L_i = L/N;
%Creación del ángulo de torsión
theta_i = zeros(1,N);
for j = 1:N
    if j < 2
        theta_i(1) = theta_1;
    else
        theta_i(j) = theta_i(j-1) + Delta_theta;
    end
end

theta_1 =     (theta_1 * pi)     / 180; %Rad
Delta_theta = (Delta_theta * pi) / 180; %Rad
theta_i =     (theta_i .* pi)    / 180; %Rad

% Densidad del material de la pala
CFRP = 1410; %kg/m^3
GFRP = 1500; %kg/m^3
GFEpoxi = 1700; %kg/m^3
dens_pala = [CFRP GFRP GFEpoxi];



%% Fórmulas para el cálculo inicial y completo de la pala de la turbina eólica.

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
% Área de cada segmento de la pala
S_i = ((c_left_i + c_right_i) / 2) * L_i;

%Definición del brazo
cateto_buje = (Buje/2) - (Punta/2);
R_brazo = sqrt(cateto_buje.^2 + L.^2);
brazo_i = (((2*i) -1) .* R_brazo) / (2*N);


% Cálculo del volumen del frustum piramidal irregular

%Me falla que estas también van variando dependiendo del segmento!!!!
ancho_buje = 2; %m
ancho_punta = 0.25; %m

area_base_menor = ancho_punta .* c_right_i;
area_base_mayor = ancho_buje .* c_left_i;
v_frustum = (L_i/3) .* (area_base_mayor + area_base_menor + sqrt(area_base_menor .* area_base_mayor));

% Masa de cada segmento de la pala
S_pala = sum(S_i); 
m_i = (S_i/S_pala) * masa_pala;

% Momento inercia del área de un trapecio
I_area = (L_i^3).*((c_right_i.^2) + (4.*c_right_i.*c_left_i) + (c_left_i.^2)) ./ (36 .* (c_right_i + c_left_i));
I_general = dens_pala(1) .* I_area;
steiner_theorem = m_i .* (brazo_i.^2);

% Momento de inercia general
I = I_general + steiner_theorem;

% Se calcula el tiempo, lo que tarda el viento para diferentes velocidades en atravesar el
% segmento
%intervalo_tiempo = c_i ./ u(v);
M = length(u_viento);
intervalo_tiempo = zeros(M,N);
for j = 1:M
    for j2 = 1:N
        intervalo_tiempo(j,j2) = c_i(j2) ./ u_viento(j);
    end
end

%% Cuando solo presenta ángulo de cabeceo
% Fuerza del viento
%F_viento_i = (1/2) .* Ro .* S_i .* (u(v).^2);

F_viento_i = zeros(M,N);
for j = 1:M
    for j2 = 1:N
        F_viento_i(j,j2) = (1/2) .* Ro .* S_i(j2) .* u_viento(j);
    end 
end
% Fuerza normal
F_normal_i = F_viento_i .* sin(theta_1);
% Momento de torsión
torque_0 = F_normal_i .* brazo_i;
% Torque global
torque_global_0 = sum(torque_0,2);



% Aceleración angular
alpha_ang_0 = torque_0 ./ I;
% Velocidad angular
Omega_0 = tiempo_analisis .* alpha_ang_0;
%Desconozco si hay que realizar la suma o la obtención del mayor
%valor de esta
Omega_0_max = max(Omega_0,[],2);
% Potencia de la pala
potencia_0 = torque_global_0 .* Omega_0_max;

%% Cuando presenta ángulo de cabeceo y torsión de los segmentos

% Fuerza normal
F_normal_i_torsion = F_viento_i .* sin(theta_i);

% Momento de torsión
%Se necesita multiplicar por el seno???? del ángulo de torsión para
%sacar el torque del área efectiva
torque_1 = zeros(M,N);
for j = 1:M
    for j2 = 1:N
        if j2 < 2
            torque_1(j,j2) = F_normal_i_torsion(j,j2) .* brazo_i(j2);
        else
            torque_1(j,j2) = F_normal_i_torsion(j,j2) .* brazo_i(j2) .* cos(Delta_theta);
        end
    end
end

% Torque global
torque_global_1 = sum(torque_1,2);
% Aceleración angular
alpha_ang_1 = torque_1 ./ I;
% Velocidad angular
Omega_1 = tiempo_analisis .* alpha_ang_1;
%Es con el max del segmento?, creo que sí.
Omega_1_max = max(Omega_1,[],2);
% Potencia de la pala
potencia_1 = torque_global_1 .* Omega_1_max;

%% Eficiencia

% Se calcula el % de mejora o empeoramiento mediante la torsión de los
% segmentos de la pala

eta = potencia_1 ./ potencia_0;

%% Representaciones

%Potencia obtenida dependiendo de la velocidad del viento
x = u_viento;
y0 = potencia_0.';
y1 = potencia_1.';
plot(x,y0);
hold on;
plot(x,y1);
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