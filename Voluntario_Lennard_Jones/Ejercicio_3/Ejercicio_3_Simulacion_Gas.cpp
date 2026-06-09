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
    int N = 169;
    double L = 32.0;

    // Genero las posiciones iniciales en una red cuadrada equiespaciada
    vector < double > x(N), y(N); // Creo los vectores para almacenar las posiciones x e y de los átomos

    int numero_lados = sqrt(N);
    double separacion = L / numero_lados;

    int i = 0;
    for (int fila = 0; fila < numero_lados; fila++)
    {
        for (int columna = 0; columna < numero_lados; columna++)
        {
            x[i] = (columna + 0.5) * separacion;
            y[i] = (fila + 0.5) * separacion;
            
            i++;
        }
    }

    // Genero un array de ángulos aleatorios entre 0 y 2pi para cada átomo
    vector < double > angulos(N);

    double pi = acos(-1.0); // Defino el valor de pi

    uniform_real_distribution<double> distribucion_angular(0.0, 2.0 * pi);

    for (int i = 0; i < N; i++) 
    {
        angulos[i] = distribucion_angular(generator);
    }

    // Genero los arrays de las velocidades en los ejes X e Y a partir de los ángulos aleatorios
    // Se deben sustituir las velocidades iniciales v = 1, 1.5, 2, 2.5, 3 para obtener todas las simulaciones

    vector < double > v_x(N), v_y(N);

    for (int i = 0; i < N; i++) 
    {
        v_x[i] = 3.0 * cos(angulos[i]);
        v_y[i] = 3.0 * sin(angulos[i]);
    }

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

    // Genero el fichero para guardar las velocidades de los átomos
    ofstream archivo_velocidades("Velocidades_Atomos_v_3_ej_3.txt");

    // Escribo las velocidades iniciales en los ficheros
    for (int i = 0; i < N; i++) 
    {
        archivo_velocidades << v_x[i] << ' ' << v_y[i] << ' ';
    }
    archivo_velocidades << '\n';

    // Inicio el ciclo for que me va a recorrer los pasos del (2) al (5)
    // Hago el ciclo for para asegurar que se completen 250 segundos

    double M = 250 / h;
    M = (int)M; // Número de pasos temporales

    // Genero el contador para almacenar el momento transferido a la caja
    double momento = 0.0;

    // Genero el array para almacenar la presión que ejercen las partículas en las paredes de la caja
    vector < double > P(M, 0.0);

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

        // Calculo el momento transferido a la caja en cada instante de tiempo por las partículas que atraviesan las fronteras
        for (int i = 0; i < N; i++) 
        {
            if (x[i] < 0.0) momento += 2 * abs(v_x[i]);
            if (x[i] >= L) momento += 2 * abs(v_x[i]);
            if (y[i] < 0.0) momento += 2 * abs(v_y[i]);
            if (y[i] >= L) momento += 2 * abs(v_y[i]);
        }

        // Aplico condiciones periódicas
        for (int i = 0; i < N; i++) {
            while (x[i] < 0.0) x[i] += L;
            while (x[i] >= L) x[i] -= L;
            while (y[i] < 0.0) y[i] += L;
            while (y[i] >= L) y[i] -= L;
        }
        // Los while se usan para los casos en los que una partícula atraviese varias veces la frontera en un solo paso temporal

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

        // Calculo la presión que ejercen las partículas en las paredes de la caja a partir de la fuerza total y el área de la caja
        P[k-1] = momento / (4 * L * t); // El área de la caja es 4 * L porque cada pared tiene una longitud L y hay 4 paredes

        // Guardo las velocidades en los ficheros
        for (int i = 0; i < N; i++) 
        {
            archivo_velocidades << v_x[i] << ' ' << v_y[i] << ' ';
        }
        archivo_velocidades << '\n';
    }

    archivo_velocidades.close();

    // Guardo los datos de la presión en un fichero
    ofstream archivo_presion("Presion_v_3_ej_3.txt");

    for (int k = 0; k < M; k++) 
    {
        archivo_presion << P[k] << '\n';
    }
    archivo_presion.close();

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