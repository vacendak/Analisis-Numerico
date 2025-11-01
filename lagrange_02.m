%% ========================================================================
%  INTERPOLACIÓN CÚBICA HACIA ADELANTE DE LAGRANGE
%  Solicita datos al usuario y genera análisis completo
%% ========================================================================

clear; close all;

%% FUNCIÓN: Interpolación de Lagrange (General)
function P = lagrange_interpolacion(x_datos, y_datos, x_eval)
    n = length(x_datos);
    P = zeros(size(x_eval));
    
    for i = 1:n
        L = ones(size(x_eval));
        for j = 1:n
            if i ~= j
                L = L .* (x_eval - x_datos(j)) / (x_datos(i) - x_datos(j));
            end
        end
        P = P + y_datos(i) * L;
    end
end

%% FUNCIÓN: Obtener expresión simbólica del polinomio
function [polinomio_str, coeficientes] = obtener_polinomio_simbolico(x_datos, y_datos)
    % Calcula la expresión explícita del polinomio de Lagrange
    syms x;
    n = length(x_datos);
    P_sim = sym(0);
    
    for i = 1:n
        L_sim = sym(1);
        for j = 1:n
            if i ~= j
                L_sim = L_sim * (x - x_datos(j)) / (x_datos(i) - x_datos(j));
            end
        end
        P_sim = P_sim + y_datos(i) * L_sim;
    end
    
    % Expandir y simplificar
    P_sim = expand(P_sim);
    P_sim = simplify(P_sim);
    
    % Convertir a string
    polinomio_str = char(P_sim);
    
    % Extraer coeficientes
    P_pol = sym2poly(P_sim);
    coeficientes = P_pol;
end

%% ========================================================================
%  PROGRAMA PRINCIPAL
%% ========================================================================

fprintf('==================================================================\n');
fprintf('  INTERPOLACIÓN CÚBICA DE LAGRANGE HACIA ADELANTE\n');
fprintf('==================================================================\n\n');

%% ENTRADA DE DATOS
fprintf('--- ENTRADA DE DATOS ---\n\n');

% Solicitar x_datos
fprintf('Introducir los valores de x (Entre corchetes y separados por espacios):\n');
x_datos = input('x_datos = ');

% Solicitar y_datos
fprintf(['\nIntroducir los valores de y (Entre corchetes y separados por ' ...
    'espacios):\n']);
y_datos = input('y_datos = ');

% Validar que tengan la misma longitud
if length(x_datos) ~= length(y_datos)
    error('Los vectores x_datos y y_datos deben tener la misma longitud');
end

% Solicitar índice de inicio
fprintf('\nÍndice del punto de partida (1 a %d):\n', length(x_datos)-3);
fprintf('Nota: Se necesitan 4 puntos consecutivos para interpolación cúbica\n');
indice_inicio = input('indice_inicio = ');

% Validar índice
if indice_inicio < 1 || indice_inicio > length(x_datos)-3
    error('Índice inválido. Debe estar entre 1 y %d', length(x_datos)-3);
end

% Solicitar punto específico a evaluar
fprintf('\nIntroducir el valor de x donde desea evaluar el polinomio:\n');
x_especifico = input('x_especifico = ');

%% PROCESAMIENTO
fprintf('==================================================================\n');
fprintf('--- PROCESAMIENTO ---\n\n');

% Seleccionar puntos para interpolación cúbica
x_usados = x_datos(indice_inicio:indice_inicio+3);
y_usados = y_datos(indice_inicio:indice_inicio+3);

% Calcular expresión simbólica del polinomio
fprintf('Calculando la expresión explícita del polinomio...\n');
[polinomio_str, coefs] = obtener_polinomio_simbolico(x_usados, y_usados);

fprintf('\n--- POLINOMIO DE INTERPOLACIÓN ---\n\n');
fprintf('P3(x) = %s\n\n', polinomio_str);

% Mostrar en forma estándar
fprintf('Forma estándar (coeficientes de mayor a menor grado):\n');
fprintf('P3(x) = ');
grado = length(coefs) - 1;
for i = 1:length(coefs)
    potencia = grado - i + 1;
    coef = coefs(i);
    
    if i == 1
        if potencia > 0
            fprintf('%.6g x^%d', coef, potencia);
        else
            fprintf('%.6g', coef);
        end
    else
        if coef >= 0
            fprintf(' + ');
        else
            fprintf(' - ');
            coef = abs(coef);
        end
        
        if potencia > 1
            fprintf('%.6g x^%d', coef, potencia);
        elseif potencia == 1
            fprintf('%.6g x', coef);
        else
            fprintf('%.6g', coef);
        end
    end
end
fprintf('\n\n');

% Evaluar en el punto específico
y_especifico = lagrange_interpolacion(x_usados, y_usados, x_especifico);

fprintf('--- EVALUACIÓN ---\n\n');
fprintf('P3(%.4f) = %.2f\n\n', x_especifico, y_especifico);

%% VISUALIZACIÓN
fprintf('--- GENERANDO GRÁFICO ---\n\n');

% Crear vector para graficar el polinomio
x_min = min(x_usados);
x_max = max(x_usados);
rango = x_max - x_min;
x_plot = linspace(x_min - 0.1*rango, x_max + 0.1*rango, 300);
y_plot = lagrange_interpolacion(x_usados, y_usados, x_plot);

% Crear figura con dimensiones especificadas
fig = figure('Position', [200, 200, 1100, 700]);

% Graficar
plot(x_datos, y_datos, 'o', 'MarkerSize', 3, 'MarkerFaceColor', [0.7 0.7 0.7], ...
    'LineWidth', 1.5, 'DisplayName', 'Datos completos');
hold on;

plot(x_usados, y_usados, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'b', ...
    'LineWidth', 2, 'DisplayName', 'Puntos usados en interpolación');

plot(x_plot, y_plot, '-', 'LineWidth', 2, 'Color', "#77AC30", ...
    'DisplayName', 'Polinomio P(x)');

plot(x_especifico, y_especifico, 'o', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', ...
    6, 'MarkerFaceColor', [0 0.8 0], 'LineWidth', 3, 'DisplayName', ...
    sprintf('Punto evaluado: (%.3f, %.2f)', x_especifico, y_especifico));

% Añadir líneas de referencia
plot([x_especifico, x_especifico], [min(ylim), y_especifico], '--','Color', ...
    [0.8500 0.3250 0.0980],'LineWidth', 1, 'HandleVisibility', 'off');
plot([min(xlim), x_especifico], [y_especifico, y_especifico], 'g--','Color', ...
    [0.8500 0.3250 0.0980],'LineWidth', 1, 'HandleVisibility', 'off');

grid on;
xlabel('N_{14C}', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Datación (años)', 'FontSize', 14, 'FontWeight', 'bold');
title({'Interpolación Cúbica de Lagrange (Hacia Adelante)', ...
       sprintf('P(x) = %s', polinomio_str)}, ...
       'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);

% Añadir cuadro de texto con información
dim = [0.15 0.72 0.25 0.15];
str = {sprintf('Puntos usados: %d', length(x_usados)), ...
       sprintf('Punto de partida: x = %.3f', x_usados(1)), ...
       sprintf('Evaluación: x = %.3f', x_especifico), ...
       sprintf('Resultado: y = %.2f', y_especifico)};
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', ...
           'BackgroundColor', 'white', 'EdgeColor', 'black', ...
           'LineWidth', 1.5, 'FontSize', 10);

hold off;

%% EXPORTAR GRÁFICO
nombre_archivo = sprintf('interpolacion_lagrange_%s.png', datetime('now', ...
    'Format','yyyyMMdd_HHmmss'));
fprintf('Exportando gráfico como: %s\n', nombre_archivo);
print(fig, nombre_archivo, '-dpng', '-r300');
fprintf('Gráfico exportado exitosamente con resolución de 300 ppp.\n\n');

%% RESUMEN FINAL
fprintf('==================================================================\n');
fprintf('--- RESUMEN ---\n');
fprintf('==================================================================\n');
fprintf('Datos de entrada: %d puntos\n', length(x_datos));
fprintf('Puntos utilizados: 4 (interpolación cúbica)\n');
fprintf('Índice de inicio: %d\n', indice_inicio);
fprintf('Rango de interpolación: [%.3f, %.3f]\n', x_usados(1), x_usados(end));
fprintf('Punto evaluado: x = %.4f\n', x_especifico);
fprintf('Resultado: y = %.2f\n', y_especifico);
fprintf('Gráfico guardado: %s\n', nombre_archivo);
fprintf('==================================================================\n');
fprintf('\nProceso completado exitosamente!\n');