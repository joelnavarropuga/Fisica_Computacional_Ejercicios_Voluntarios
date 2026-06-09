# include <iostream>
# include <cmath>
# include <vector>
# include <fstream>
# include <chrono>
# include <random>
# include <algorithm>
# include <string>

using namespace std;

void aceleracion(const vector < double > &x, const vector < double > &y, int N, double L, vector < double > &a_x, vector < double > &a_y);

int main()
{
    auto inicio = chrono::high_resolution_clock::now(); // Inicio del cronómetro para medir el tiempo de ejecución

    unsigned seed = chrono::system_clock::now().time_since_epoch().count();  // Inicializo la semilla para generar números aleatorios
    mt19937_64 generator(seed);  // Inicializo el generador de números aleatorios

    // Defino el número de átomos y el tamaño de la caja
    int N = 144;
    double L = 12.0;

    // Genero las posiciones iniciales en una red cuadrada equiespaciada
    vector < double > x(N), y(N); // Creo los vectores para almacenar las posiciones x e y de los átomos

    int numero_lados = sqrt(N);

    int i = 0;
    for (int fila = 0; fila < numero_lados; fila++)
    {
        for (int columna = 0; columna < numero_lados; columna++)
        {
            x[i] = columna;
            y[i] = fila;
            
            i++;
        }
    }

    // Genero los arrays de las velocidades nulas

    vector < double > v_x(N, 0.0), v_y(N, 0.0);

    /////////////////////////////////////////////////////////
    /// ALGORITMO DE VERLET ///
    /////////////////////////////////////////////////////////

    // Inicializo t a 0
    double t = 0.0;

    // Defino el paso temporal h
    double h = 0.002;

    // Calculo la aceleración inicial de cada átomo a partir de las posiciones iniciales
    vector < double > a_x(N, 0.0), a_y(N, 0.0);
    aceleracion(x, y, N, L, a_x, a_y);

    // Genero los ficheros donde voy a guardar las posiciones y velocidades de los átomos
    // Las posiciones se van a guardar en forma x1 y1 x2 y2 ... xN yN
    // Las velocidades se van a guardar en forma vx1 vy1 vx2 vy2 ... vxN vyN

    ofstream archivo_posiciones("Posiciones_Atomos_ej_7.txt");
    ofstream archivo_velocidades("Velocidades_Atomos_ej_7.txt");
    ofstream archivo_posicion_atomo("Posicion_Par_Atomos_ej_7.txt");

    // Escribo las posiciones y velocidades iniciales en los ficheros
    for (int i = 0; i < N; i++) 
    {
        archivo_posiciones << x[i] << ' ' << y[i] << ' ';
        archivo_velocidades << v_x[i] << ' ' << v_y[i] << ' ';
    }
    archivo_posiciones << '\n';
    archivo_velocidades << '\n';

    // Genero dos variables para guardar las posiciones del átomo 0 y 1, sin aplicar las condiciones periódicas en el fichero
    double x_0_aux = x[0];
    double y_0_aux = y[0];
    double x_1_aux = x[1];
    double y_1_aux = y[1];

    // Escribo las posiciones de los átomos 0 y 1en el fichero
    archivo_posicion_atomo << x_0_aux << ' ' << y_0_aux << ' ' << x_1_aux << ' ' << y_1_aux << '\n';

    // Inicio el ciclo for que me va a recorrer los pasos del (2) al (5)
    // Hago el ciclo for para asegurar que se completen 250 segundos

    double M = 3500 / h;
    M = (int)M; // Número de pasos temporales

    // Inicio el ciclo temporal
    for (int k = 1; k <= M; k++) 
    {
        /// PASO 2. Calculo r_i (t + h) y w_i

        vector < double > w_x(N), w_y(N);

        for (int i = 0; i < N; i++) 
        {
            w_x[i] = v_x[i] + (h / 2.0) * a_x[i];
            w_y[i] = v_y[i] + (h / 2.0) * a_y[i];
        }

        // Calculo la nueva posición en t + h
        for (int i = 0; i < N; i++) 
        {
            x[i] = x[i] + h * w_x[i];
            y[i] = y[i] + h * w_y[i];
        }

        // Trayectoria real de los átomos 0 y 1 sin aplicar condiciones periódicas
        x_0_aux += h * w_x[0];
        y_0_aux += h * w_y[0];
        x_1_aux += h * w_x[1];
        y_1_aux += h * w_y[1];

        archivo_posicion_atomo << x_0_aux << ' ' << y_0_aux << ' ' << x_1_aux << ' ' << y_1_aux << '\n';

        // Aplico condiciones periódicas
        for (int i = 0; i < N; i++) {
            if (x[i] < 0.0) x[i] += L;
            if (x[i] >= L) x[i] -= L;
            if (y[i] < 0.0) y[i] += L;
            if (y[i] >= L) y[i] -= L;
        }

        /// PASO 3. Calculo las aceleraciones en t + h usando las posiciones en t + h

        for (int i = 0; i < N; i++) {
            a_x[i] = 0.0;
            a_y[i] = 0.0;
        }

        aceleracion(x, y, N, L, a_x, a_y);

        /// PASO 4. Evalúo las velocidades en t + h

        for (int i = 0; i < N; i++) 
        {
            v_x[i] = w_x[i] + (h / 2) * a_x[i];
            v_y[i] = w_y[i] + (h / 2) * a_y[i];
        }

        /// PASO 5. Poner t = t + h

        t += h;

        // Guardo las posiciones y velocidades en los ficheros

        for (int i = 0; i < N; i++) 
        {
            archivo_posiciones << x[i] << ' ' << y[i] << ' ';
            archivo_velocidades << v_x[i] << ' ' << v_y[i] << ' ';
        }
        archivo_posiciones << '\n';
        archivo_velocidades << '\n';

        // Reescalo las velocidades para aumentar la temperatura en t = 500, 1000, 2000, 3000 segundos
        if ((k == 250000) || (k == 500000) || (k == 1000000) || (k == 1500000))
        {
            for (int i = 0; i < N; i++) 
            {
                v_x[i] *= 1.1;
                v_y[i] *= 1.1;
            }
        }
    }

    archivo_posiciones.close();
    archivo_velocidades.close();
    archivo_posicion_atomo.close();
    
    auto fin = chrono::high_resolution_clock::now(); // Fin del cronómetro para medir el tiempo de ejecución
    double tiempo_seg = chrono::duration<double>(fin - inicio).count();
    
    cout << "Tiempo de ejecucion: " << tiempo_seg << " s" << endl;

    return 0;
}

// Defino la función que calcula la aceleración inicial de cada átomo
void aceleracion(const vector < double > &x, const vector < double > &y, int N, double L, vector < double > &a_x, vector < double > &a_y)  // const es para indicar que los vectores x e y no se van a modificar dentro de la función, lo que optimiza el código al evitar copias innecesarias de los vectores
{
    for (int i = 0; i < N; i++) 
    {
        for (int j = i + 1; j < N; j++) // con j = i + 1 evito calcular dos veces la interacción entre los átomos i y j
        {
            // Calculo la distancia entre los átomos i y j teniendo en cuenta las condiciones periódicas de la caja
            double x_ij = x[i] - x[j];
            if (x_ij > L / 2.0) x_ij -= L;
            if (x_ij < -L / 2.0) x_ij += L;

            double y_ij = y[i] - y[j];
            if (y_ij > L / 2.0) y_ij -= L;
            if (y_ij < -L / 2.0) y_ij += L;

            double r_ij = sqrt(x_ij * x_ij + y_ij * y_ij);

            if ((r_ij > 0.0) && (r_ij < 3.0)) // No hay interacción entre los átomos si la distancia es mayor que 3.0
            {
                double inv_r_ij = 1.0 / r_ij;
                double inv_r_ij_8 = pow(inv_r_ij, 8);
                double inv_r_ij_14 = pow(inv_r_ij, 14);

                double f = 24.0 * (2.0 * inv_r_ij_14 - inv_r_ij_8);
                double fx = f * x_ij;
                double fy = f * y_ij;

                a_x[i] += fx;
                a_y[i] += fy;

                // Por la tercera ley de Newton, la fuerza que el átomo j ejerce sobre el átomo i es igual y opuesta a la fuerza que el átomo i ejerce sobre el átomo j
                a_x[j] -= fx;
                a_y[j] -= fy;
            }
        }
    }
}