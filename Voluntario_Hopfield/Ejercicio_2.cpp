# include <iostream>
# include <math.h>
# include <vector>
# include <fstream>
# include <chrono>
# include <random>
# include <algorithm>
# include <string>

using namespace std;

vector < vector < int > > configuracion_inicial_aleatoria(int N, mt19937_64& generator);
vector < vector < int > > configuracion_inicial_deformada(vector < vector < int > > patron_original, int N, mt19937_64& generator);
vector < vector < vector < vector < double > > > > omega(vector < vector < int > > patron_original, int N);
void algoritmo_Metropolis(vector < vector < int > > s, vector < vector < vector < vector < double > > > > omega_matriz, int N, mt19937_64& generator, ofstream& fichero_neuronas, double T);

int main()
{
    auto inicio = chrono::high_resolution_clock::now(); // Inicio del cronómetro para medir el tiempo de ejecución

    unsigned seed = chrono::system_clock::now().time_since_epoch().count();  // Inicializo la semilla para generar números aleatorios
    mt19937_64 generator(seed);  // Inicializo el generador de números aleatorios

    // Leo el patrón de entrada desde un archivo
    ifstream archivo_lectura("Patron_original.txt");

    // Creo una matriz para almacenar el patrón de entrada
    int N = 150;
    vector < vector < int > > patron_original(N, vector < int > (N, 0)); 

    // Leo el patrón de entrada desde el archivo y lo almaceno en la matriz
    ifstream archivo_escritura("Patron_original.txt");
    
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++) 
        {
            archivo_escritura >> patron_original[i][j];
        }
    }

    archivo_escritura.close(); // Cierro el archivo

    // Creo la matriz omega de interacción entre las neuronas llamando a la función omega
    vector < vector < vector < vector < double > > > > omega_matriz = omega(patron_original, N);

    // Genero una configuración inicial aleatoria para la red de neuronas, con valores de 0 o 1 al azar, llamando a la función
    vector < vector < int > > s_aleatorio = configuracion_inicial_aleatoria(N, generator);

    // Creo un array con las temperaturas a las que voy a ejecutar el algoritmo de Metropolis
    vector < double > T = {0.001, 0.01, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1};
    
    // Ejecuto el algortimo de Metropolis para una configuración inicial aleatoria
    for (int i = 0; i < T.size(); i++)
    {
        ofstream fichero_neuronas_aleatorias("Configuracion_inicial_neuronas_aleatoria_T_" + to_string(T[i]) + ".txt");
    
        algoritmo_Metropolis(s_aleatorio, omega_matriz, N, generator, fichero_neuronas_aleatorias, T[i]);
    }

    // Genero una configuración inicial parecida al patrón de entrada
    vector < vector < int > > s_deformado = configuracion_inicial_deformada(patron_original, N, generator);

    // Ejecuto el algoritmo de Metropolis para una configuración inicial parecida al patrón de entrada
    for (int i = 0; i < T.size(); i++)
    {
        ofstream fichero_neuronas_deformado("Configuracion_inicial_neuronas_deformada_T_" + to_string(T[i]) + ".txt");

        algoritmo_Metropolis(s_deformado, omega_matriz, N, generator, fichero_neuronas_deformado, T[i]);
    }

    auto fin = chrono::high_resolution_clock::now(); // Fin del cronómetro para medir el tiempo de ejecución
    
    double tiempo_seg = chrono::duration<double>(fin - inicio).count();
    cout << "Tiempo de ejecucion: " << tiempo_seg << " s" << endl;

    return 0;
}

// Función para generar una configuración inicial aleatoria para la red de neuronas, con valores de 0 o 1 al azar
vector < vector < int > > configuracion_inicial_aleatoria(int N, mt19937_64& generator)
{
    uniform_int_distribution < int > bit_dist(0, 1);
    vector < vector < int > > s(N, vector < int > (N, 0));

    for (int i = 0; i < N; i++) 
    {
        for (int j = 0; j < N; j++) 
        {
            s[i][j] = bit_dist(generator);  // Asigno un valor de 0 o 1 al azar a cada neurona
        }
    }
    return s;
}

vector < vector < int > > configuracion_inicial_deformada(vector < vector < int > > patron_original,int N, mt19937_64& generator)
{    
    vector < vector < int > > s(N, vector < int > (N, 0));
    s = patron_original;

    // Deformo el patrón inicial cambiando el valor de activación de un 10% de las neuronas al azar

    uniform_int_distribution < int > distribution(0, N - 1);
    
    for (int k = 0; k < N * N / 10; k++)
    {
        int i = distribution(generator);
        int j = distribution(generator);

        if (s[i][j] == 0) 
        {
            s[i][j] = 1;  // Cambio el valor de activación de la neurona seleccionada a 1
        } 

        else if (s[i][j] == 1)
        {
            s[i][j] = 0;  // Cambio el valor de activación de la neurona seleccionada a 0
        }
    }

    return s;
}

// Función para rellenar la matriz omega de 4 dimensiones 
vector < vector < vector < vector < double > > > > omega(vector < vector < int > > patron_original, int N)
{
    double sum = 0.0;
    double a;

    // Defino la a^mu correspondiente al patrón de entrada
    for (int r = 0; r < N; r++) 
    {
        for (int t = 0; t < N; t++) 
        {
           sum += patron_original[r][t];
        }
    }

    a = sum / (N * N);

    // Construyo omega de interacción entre la neurona i,j y la neurona k,l 
    vector < vector < vector < vector < double > > > > omega_matriz(N, vector < vector < vector < double > > >(N, vector < vector < double > >(N, vector < double >(N, 0.0))));

    for (int i = 0; i < N; i++) 
    {
        for (int j = 0; j < N; j++) 
        {
            for (int k = 0; k < N; k++) 
            {
                for (int l = 0; l < N; l++) 
                {
                    if ((i == k) && (j == l)) 
                    {
                        omega_matriz[i][j][k][l] = 0.0;  // No hay auto-conexiones
                    } 

                    else 
                    {
                        omega_matriz[i][j][k][l] = (1.0 / (N * N)) * (patron_original[i][j] - a) * (patron_original[k][l] - a);
                    }
                }
            }
        }
    }
    return omega_matriz;
}

// Función para aplicar el algoritmo de Metropolis
void algoritmo_Metropolis(vector < vector < int > > s, vector < vector < vector < vector < double > > > > omega_matriz, int N, mt19937_64& generator, ofstream& fichero_neuronas, double T)
{
    //PASO 0: Asigno un valor inicial a la temperatura y creo una configuración inicial de la red de neuronas

    // Guardo la configuración inicial de las neuronas
    for (int i = 0; i < N; i++) 
    {
        for (int j = 0; j < N; j++) 
        {
            fichero_neuronas << s[i][j] << " ";
        }
        fichero_neuronas << endl;
    }
    fichero_neuronas << endl;

    // Defino un paso Monte Carlo
    int paso_Monte_Carlo = N * N;

    // Defino el numero de pasos Monte Carlo que voy a realizar
    int num_pasos = 20 * paso_Monte_Carlo;

    for (int k = 0; k < num_pasos; k++) 
    {

        //PASO 1: Elijo un punto al azar de la red 

        uniform_int_distribution < int > distribution(0, N - 1);

        int i = distribution(generator);
        int j = distribution(generator);

        //PASO 2: Calculo la probabilidad de transición del estado inicial a un estado con el valor de activación (0, 1) de la neurona seleccionada cambiado """

        // Calculo la variación de energía asociada a la energía de interacción
        double dE_interaccion = 0.0;
        double theta = 0.0;

        if (s[i][j] == 0) 
        {
            for (int k = 0; k < N; k++) 
            {
                for (int l = 0; l < N; l++) 
                {
                    dE_interaccion -= 2 * omega_matriz[i][j][k][l] * s[k][l];
                    theta += omega_matriz[i][j][k][l];
                }
            }
        }

        else if (s[i][j] == 1)
        {
            for (int k = 0; k < N; k++) 
            {
                for (int l = 0; l < N; l++) 
                {
                    dE_interaccion += 2 * omega_matriz[i][j][k][l] * s[k][l];
                    theta += omega_matriz[i][j][k][l];
                }
            }
        }

        theta = 0.5 * theta;

        // Calculo la variación de energía asociada a la theta
        double dE_theta;

        if (s[i][j] == 0) 
        {
            dE_theta = theta;
        } 

        else if (s[i][j] == 1)
        {
            dE_theta = -theta;
        }

        // Calculo la variación total de energía 
        double dE;

        dE = dE_interaccion + dE_theta;

        // Defino la probabilidad de transición
        double P;

        P = min(1.0, exp(-dE / T));

        //PASO 3: Genero un número aleatorio entre 0 y 1 y comparo con la probabilidad de transición para aceptar o rechazar el cambio de estado """

        //Genero un número aleatorio entre 0 y 1
        uniform_real_distribution < double > r_distribution(0.0, 1.0);

        double epsilon = r_distribution(generator);

        if (epsilon < P)
        {
            s[i][j] = 1 - s[i][j];  // Cambio el valor de activación de la neurona seleccionada al valor opuesto (0 a 1 o 1 a 0)
        }

        // Guardo la configuración final de las neuronas cada paso Monte Carlo
        if ((k + 1) % paso_Monte_Carlo == 0) 
        {
            for (int i = 0; i < N; i++) 
            {
                for (int j = 0; j < N; j++) 
                {
                    fichero_neuronas << s[i][j] << " ";
                }
                fichero_neuronas << endl;
            }
            fichero_neuronas << endl;
        }
    }

    fichero_neuronas.close(); // Cierro el fichero donde he guardado la configuración de las neuronas
}