% ====================================================================
% METODO DE INTERPOLACION POR DIFERENCIAS PROGRESIVAS DE NEWTON
% ====================================================================
% Este programa calcula el polinomio interpolador usando el metodo de
% diferencias progresivas (hacia adelante) de Newton para puntos
% igualmente espaciados.
%
% ENTRADA:
%   - Vector de valores x (puntos igualmente espaciados)
%   - Vector de valores y (valores de la funcion en esos puntos)
%   - Valor del punto x donde se desea interpolar
%
% SALIDA:
%   - Tabla de diferencias progresivas
%   - Expresion simbolica del polinomio interpolador
%   - Valor interpolado en el punto solicitado
% ====================================================================

clear;
format long;

% --------------------------------------------------------------------
% SECCION 1: ENTRADA DE DATOS
% --------------------------------------------------------------------
disp('========================================================');
disp('  INTERPOLACION POR DIFERENCIAS PROGRESIVAS DE NEWTON  ');
disp('========================================================');
disp(' ');

% Solicitar los puntos x_i (deben estar igualmente espaciados)
disp('Introduce el vector de valores x_i (entre corchetes, separados por espacios):');
disp('Ejemplo: [0 1 2 3 4]');
vector_x = input('x_i = ');

% Solicitar los valores y_i correspondientes
disp(' ');
disp('Introduce el vector de valores y_i (entre corchetes, separados por espacios):');
disp('Ejemplo: [1 3 9 19 33]');
vector_y = input('y_i = ');

% Solicitar el punto donde interpolar
disp(' ');
punto_interpolar = input('Introduce el valor de x donde deseas interpolar: ');

disp(' ');
disp('--------------------------------------------------------');

% --------------------------------------------------------------------
% SECCION 2: VALIDACION DE DATOS
% --------------------------------------------------------------------

% Verificar que los vectores tengan la misma longitud
if length(vector_x) ~= length(vector_y)
    error('ERROR: Los vectores x e y deben tener la misma longitud.');
end

% Verificar que haya al menos 2 puntos
n_puntos = length(vector_x);
if n_puntos < 2
    error('ERROR: Se necesitan al menos 2 puntos para interpolar.');
end

% Verificar que los puntos esten igualmente espaciados
espaciado = diff(vector_x);
tolerancia = 1e-10;
if max(abs(espaciado - espaciado(1))) > tolerancia
    error('ERROR: Los puntos x_i deben estar igualmente espaciados.');
end

% Calcular el espaciado h
h = espaciado(1);
fprintf('Espaciado constante detectado: h = %.6f\n', h);
disp(' ');

% --------------------------------------------------------------------
% SECCION 3: CONSTRUCCION DE LA TABLA DE DIFERENCIAS PROGRESIVAS
% --------------------------------------------------------------------

% Inicializar la matriz de diferencias
% Columnas: x_i, y_i, Delta_y, Delta2_y, Delta3_y, ...
tabla_diferencias = zeros(n_puntos, n_puntos + 1);

% Primera columna: valores de x
tabla_diferencias(:, 1) = vector_x';

% Segunda columna: valores de y (diferencias de orden 0)
tabla_diferencias(:, 2) = vector_y';

% Calcular las diferencias progresivas de orden superior
for orden = 1:(n_puntos - 1)
    for fila = 1:(n_puntos - orden)
        % Delta^k y_i = Delta^(k-1) y_(i+1) - Delta^(k-1) y_i
        tabla_diferencias(fila, orden + 2) = ...
            tabla_diferencias(fila + 1, orden + 1) - ...
            tabla_diferencias(fila, orden + 1);
    end
end

% --------------------------------------------------------------------
% SECCION 4: MOSTRAR LA TABLA DE DIFERENCIAS
% --------------------------------------------------------------------

disp('TABLA DE DIFERENCIAS PROGRESIVAS:');
disp('--------------------------------------------------------');

% Crear encabezados de la tabla
encabezado = '   x_i        y_i    ';
for i = 1:(n_puntos - 1)
    if i == 1
        encabezado = [encabezado, sprintf('   Delta y   ')];
    else
        encabezado = [encabezado, sprintf(' Delta^%d y  ', i)];
    end
end
disp(encabezado);
disp('--------------------------------------------------------');

% Mostrar cada fila de la tabla
for fila = 1:n_puntos
    linea = sprintf('%8.4f  %10.6f', tabla_diferencias(fila, 1), ...
                    tabla_diferencias(fila, 2));
    for col = 3:(n_puntos + 1)
        if fila <= (n_puntos - col + 2)
            linea = [linea, sprintf('  %10.6f', tabla_diferencias(fila, col))];
        else
            linea = [linea, sprintf('      -      ')];
        end
    end
    disp(linea);
end
disp('--------------------------------------------------------');
disp(' ');

% --------------------------------------------------------------------
% SECCION 5: CONSTRUCCION DEL POLINOMIO INTERPOLADOR
% --------------------------------------------------------------------

% Definir la variable simbolica
syms s;

% Inicializar el polinomio con y_0
polinomio_s = tabla_diferencias(1, 2);

% Construir el polinomio usando la formula de Newton
% P(s) = y_0 + C(s,1)*Delta_y_0 + C(s,2)*Delta2_y_0 + ...
% donde C(s,k) es el coeficiente binomial generalizado

for k = 1:(n_puntos - 1)
    % Calcular el coeficiente binomial C(s,k) = s(s-1)(s-2)...(s-k+1)/k!
    coef_binomial = 1;
    for j = 0:(k - 1)
        coef_binomial = coef_binomial * (s - j);
    end
    coef_binomial = coef_binomial / factorial(k);
    
    % Agregar el termino al polinomio
    diferencia_k = tabla_diferencias(1, k + 2);
    polinomio_s = polinomio_s + coef_binomial * diferencia_k;
end

% Simplificar el polinomio en terminos de s
polinomio_s = expand(polinomio_s);

% Convertir de s a x usando s = (x - x_0) / h
syms x;
x_0 = vector_x(1);
s_en_x = (x - x_0) / h;
polinomio_x = subs(polinomio_s, s, s_en_x);
polinomio_x = expand(polinomio_x);
polinomio_x = simplify(polinomio_x);

% --------------------------------------------------------------------
% SECCION 6: MOSTRAR EL POLINOMIO INTERPOLADOR
% --------------------------------------------------------------------

disp('POLINOMIO INTERPOLADOR:');
disp('--------------------------------------------------------');
fprintf('En terminos de s = (x - %.4f) / %.4f :\n', x_0, h);
fprintf('P(s) = %s\n', char(polinomio_s));
disp(' ');
fprintf('En terminos de x:\n');
fprintf('P(x) = %s\n', char(polinomio_x));
disp('--------------------------------------------------------');
disp(' ');

% --------------------------------------------------------------------
% SECCION 7: INTERPOLACION EN EL PUNTO SOLICITADO
% --------------------------------------------------------------------

% Calcular s para el punto de interpolacion
s_interpolar = (punto_interpolar - x_0) / h;

% Evaluar el polinomio en el punto solicitado
valor_interpolado = double(subs(polinomio_x, x, punto_interpolar));

% Mostrar resultado
disp('RESULTADO DE LA INTERPOLACION:');
disp('--------------------------------------------------------');
fprintf('Punto a interpolar: x = %.6f\n', punto_interpolar);
fprintf('Valor de s: s = %.6f\n', s_interpolar);
fprintf('Valor interpolado: y(%.6f) = %.10f\n', punto_interpolar, ...
        valor_interpolado);
disp('--------------------------------------------------------');
disp(' ');

% --------------------------------------------------------------------
% SECCION 8: VERIFICACION (Opcional)
% --------------------------------------------------------------------

% Verificar que el polinomio pasa por todos los puntos originales
disp('VERIFICACION (el polinomio debe pasar por todos los puntos):');
disp('--------------------------------------------------------');
disp('   x_i      y_i dado    P(x_i)      Error');
disp('--------------------------------------------------------');

max_error = 0;
for i = 1:n_puntos
    valor_calculado = double(subs(polinomio_x, x, vector_x(i)));
    error_abs = abs(valor_calculado - vector_y(i));
    max_error = max(max_error, error_abs);
    fprintf('%8.4f  %10.6f  %10.6f  %e\n', vector_x(i), vector_y(i), ...
            valor_calculado, error_abs);
end
disp('--------------------------------------------------------');
fprintf('Error maximo: %e\n', max_error);
disp(' ');

% --------------------------------------------------------------------
% SECCION 9: GRAFICA (Opcional)
% --------------------------------------------------------------------

disp('Generando grafica...');

% Generar puntos para graficar el polinomio
x_min = min(vector_x);
x_max = max(vector_x);
rango = x_max - x_min;
x_grafica = linspace(x_min - 0.1*rango, x_max + 0.1*rango, 500);
y_grafica = double(subs(polinomio_x, x, x_grafica));

% Crear figura con dimensiones especificadas
fig = figure('Position', [200, 200, 1100, 700]);

% Graficar el polinomio
plot(x_grafica, y_grafica, '-', 'LineWidth', 2, 'Color', "#77AC30", ...
    'DisplayName', 'Polinomio P(x)');
hold on;

% Graficar los puntos originales
plot(vector_x, vector_y, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'b', ...
    'LineWidth', 2, 'DisplayName', 'Puntos usados en interpolación');

% Graficar el punto interpolado
plot(punto_interpolar, valor_interpolado, 'o', 'Color', [0.8500 0.3250 0.0980], ...
    'MarkerSize', ...
    6, 'MarkerFaceColor', [0 0.8 0], 'LineWidth', 3, 'DisplayName', ...
    sprintf('Punto evaluado: (%.6f, %.6f)', punto_interpolar, punto_interpolar));

% Añadir líneas de referencia
plot([punto_interpolar, punto_interpolar], [min(ylim), valor_interpolado], ...
    '--','Color', [0.8500 0.3250 0.0980],'LineWidth', 1, 'HandleVisibility', 'off');
plot([min(xlim), punto_interpolar], [valor_interpolado, valor_interpolado], ...
    'g--','Color', [0.8500 0.3250 0.0980],'LineWidth', 1, 'HandleVisibility', 'off');

% Configurar la grafica
grid on;
xlabel('N_{14C}', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Datación (años)', 'FontSize', 14, 'FontWeight', 'bold');
title({'Interpolacion por Diferencias Progresivas de Newton', ...
       sprintf('P(x) = %s', polinomio_x)}, ...
       'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);

% Añadir cuadro de texto con información
dim = [0.15 0.72 0.25 0.15];
str = {sprintf('Puntos usados: %d', length(vector_x)), ...
       sprintf('Evaluación: x = %.6f', punto_interpolar), ...
       sprintf('Resultado: y = %.6f', valor_interpolado)};
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', ...
           'BackgroundColor', 'white', 'EdgeColor', 'black', ...
           'LineWidth', 1.5, 'FontSize', 10);

hold off;

%% EXPORTAR GRÁFICO
nombre_archivo = sprintf('dif_prog_newton_%s.png', datetime('now', ...
    'Format','yyyyMMdd_HHmmss'));
fprintf('Exportando gráfico como: %s\n', nombre_archivo);
print(fig, nombre_archivo, '-dpng', '-r300');
fprintf('Gráfico exportado exitosamente con resolución de 300 ppp.\n\n');

disp('Grafica generada exitosamente.');
disp(' ');
disp('========================================================');
disp('              PROCESO COMPLETADO                       ');
disp('========================================================');