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
void almacenar_histograma(const vector < double > &x, const vector < double > &y, int N, double L, double r_min, double r_max, double ancho_bin, vector < int > &histograma);

int main()
{
    auto inicio = chrono::high_resolution_clock::now(); // Inicio del cronómetro para medir el tiempo de ejecución

    unsigned seed = chrono::system_clock::now().time_since_epoch().count();  // Inicializo la semilla para generar números aleatorios
    mt19937_64 generator(seed);  // Inicializo el generador de números aleatorios

    // Defino el número de átomos y el tamaño de la caja
    int N = 144;
    double L = 24.0; // Aumento el lado de la caja para reducir la densidad y hacerlo más diluido

    // Defino el histograma 
    int num_bins = 40;
    double r_min = 1.0;
    double r_max = 3.0;
    double ancho_bin = (r_max - r_min) / num_bins;
    vector < int > histograma_gas(num_bins, 0); // Inicializo el histograma con ceros

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

    // Genero un array de ángulos aleatorios entre 0 y 2pi para cada átomo
    vector < double > angulos(N);

    double pi = acos(-1.0); // Defino el valor de pi

    uniform_real_distribution<double> distribucion_angular(0.0, 2.0 * pi);

    for (int i = 0; i < N; i++) 
    {
        angulos[i] = distribucion_angular(generator);
    }
    
    // Genero los arrays de las velocidades en los ejes X e Y a partir de los ángulos aleatorios
    // El módulo de la velocidad se pone a 3 para aumentar la temperatura y que el sistema esté en estado gaseoso

    vector < double > v_x(N, 0.0), v_y(N, 0.0);

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
    double h = 0.001;

    // Calculo la aceleración inicial de cada átomo a partir de las posiciones iniciales
    vector < double > a_x(N, 0.0), a_y(N, 0.0);
    aceleracion(x, y, N, L, a_x, a_y);

    // Inicio el ciclo for que me va a recorrer los pasos del (2) al (5)
    // Hago el ciclo for para asegurar que se completen 250 segundos

    double M = 500 / h;
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

        // Almaceno en el histograma las partículas que están a una distancia entre 1 y 3 de la partícula 0 cada 100 pasos temporales 

        // Guardo en el histograma de gas partir de t = 100, para asegurar que el sistema esté en equilibrio
        if ((k >= 100000) && (k % 100 == 0))
        {
            almacenar_histograma(x, y, N, L, r_min, r_max, ancho_bin, histograma_gas);
        }
    }

    ofstream archivo_histograma_gas("Histograma_gas_ej_8.txt");
    
    for (int i = 0; i < num_bins; i++) 
    {
        archivo_histograma_gas << histograma_gas[i] << '\n';
    }
    archivo_histograma_gas.close();

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

// Función para almacenar en el histograma las partículas
void almacenar_histograma(const vector < double > &x, const vector < double > &y, int N, double L, double r_min, double r_max, double ancho_bin, vector < int > &histograma)
{
    int i = 0; // Tomamos la partícula 0 como origen
    for (int j = 1; j < N; j++)
    {
        // Calculo la distancia entre la partícula i y la j teniendo en cuenta las condiciones periódicas de la caja
        double x_ij = x[i] - x[j];
        if (x_ij > L / 2.0) x_ij -= L;
        if (x_ij < -L / 2.0) x_ij += L;

        double y_ij = y[i] - y[j];
        if (y_ij > L / 2.0) y_ij -= L;
        if (y_ij < -L / 2.0) y_ij += L;

        double r_ij = sqrt(x_ij * x_ij + y_ij * y_ij);

        // Si la distancia está en el rango 1-3, la almaceno en el bin
        if ((r_ij >= r_min) && (r_ij < r_max)) 
        {
            int bin = (int)((r_ij - r_min) / ancho_bin);

            histograma[bin]++;
        }
    }
}