# include <iostream>
# include <math.h>
# include <vector>
# include <fstream>
# include <chrono>
# include <random>
# include <algorithm>
# include <string>

using namespace std;

using namespace std;

vector < vector < int > > configuracion_inicial_aleatoria(int N, mt19937_64& generator);
vector < vector < int > > configuracion_inicial_deformada(vector < vector < vector < int > > > patron_original_mu, int N, mt19937_64& generator);
vector < vector < vector < vector < double > > > > omega(vector < vector < vector < int > > > patron_original_mu, int N, int num_patrones);
void algoritmo_Metropolis(vector < vector < int > > s, vector < vector < vector < vector < double > > > > omega_matriz, int N, mt19937_64& generator, ofstream& fichero_neuronas);

int main()
{
    auto inicio = chrono::high_resolution_clock::now(); // Inicio del cronómetro para medir el tiempo de ejecución

    unsigned seed = chrono::system_clock::now().time_since_epoch().count();  // Inicializo la semilla para generar números aleatorios
    mt19937_64 generator(seed);  // Inicializo el generador de números aleatorios

    // Leo los patrones de entrada desde los archivos
    ifstream archivo_0("0.txt");
    ifstream archivo_1("1.txt");
    ifstream archivo_2("2.txt");
    ifstream archivo_3("3.txt");
    ifstream archivo_4("4.txt");
    ifstream archivo_5("5.txt");
    ifstream archivo_6("6.txt");
    ifstream archivo_7("7.txt");
    ifstream archivo_8("8.txt");
    ifstream archivo_9("9.txt");

    // Creo las matrices para almacenar el patrón de entrada
    int N = 30;
    int num_patrones = 3;

    vector < vector < vector < int > > > patron_original_mu(num_patrones, vector < vector < int > > (N, vector < int > (N, 0)));

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++) 
        {
            archivo_0 >> patron_original_mu[0][i][j];
            archivo_1 >> patron_original_mu[1][i][j];
            archivo_2 >> patron_original_mu[2][i][j];
            //archivo_3 >> patron_original_mu[3][i][j];
            //archivo_4 >> patron_original_mu[4][i][j];
            //archivo_5 >> patron_original_mu[5][i][j];
            //archivo_6 >> patron_original_mu[6][i][j];
            //archivo_7 >> patron_original_mu[7][i][j];
            //archivo_8 >> patron_original_mu[8][i][j];
            //archivo_9 >> patron_original_mu[9][i][j];
        }
    }

    archivo_0.close();
    archivo_1.close();
    archivo_2.close();
    archivo_3.close();
    archivo_4.close();
    archivo_5.close();
    archivo_6.close();
    archivo_7.close();
    archivo_8.close();
    archivo_9.close();

    // Creo la matriz omega de interacción entre las neuronas llamando a la función omega
    vector < vector < vector < vector < double > > > > omega_matriz = omega(patron_original_mu, N, num_patrones);

    // Genero una configuración inicial aleatoria para la red de neuronas, con valores de 0 o 1 al azar, llamando a la función
    vector < vector < int > > s_aleatorio = configuracion_inicial_aleatoria(N, generator);

    // Genero el archivo donde voy a guardar la evolución de las neuronas con configuración inicial aleatoria
    ofstream fichero_neuronas_aleatorias("Configuracion_inicial_neuronas_aleatoria_mu.txt");
    
    // Ejecuto el algortimo de Metropolis para una configuración inicial aleatoria
    algoritmo_Metropolis(s_aleatorio, omega_matriz, N, generator, fichero_neuronas_aleatorias);

    // Genero una configuración inicial parecida al patrón de entrada
    vector < vector < int > > s_deformado = configuracion_inicial_deformada(patron_original_mu, N, generator);

    // Genero el archivo donde voy a guardar la configuración de las neuronas con configuración inicial parecida al patrón de entrada
    ofstream fichero_neuronas_deformado("Configuracion_inicial_neuronas_deformada_mu.txt");

    // Ejecuto el algoritmo de Metropolis para una configuración inicial parecida al patrón de entrada
    algoritmo_Metropolis(s_deformado, omega_matriz, N, generator, fichero_neuronas_deformado);

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

vector < vector < int > > configuracion_inicial_deformada(vector < vector < vector < int > > > patron_original_mu, int N, mt19937_64& generator)
{    
    vector < vector < int > > s(N, vector < int > (N, 0));
    s = patron_original_mu[0];  // Inicializo s con el primer patrón de entrada

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
vector < vector < vector < vector < double > > > > omega(vector < vector < vector < int > > > patron_original_mu, int N, int num_patrones)
{
    vector < double >  a(num_patrones, 0.0);

    // Defino la a^mu correspondiente al patrón de entrada
    for (int mu = 0; mu < num_patrones; mu++) 
    {
        double sum = 0.0;
        for (int i = 0; i < N; i++) 
        {
            for (int j = 0; j < N; j++) 
            {
                sum += patron_original_mu[mu][i][j];
            }
        }

        a[mu] = sum / (N * N);
    }

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
                        double sum_omega = 0.0;
                        for (int mu = 0; mu < num_patrones; mu++) 
                        {
                            sum_omega += (patron_original_mu[mu][i][j] - a[mu]) * (patron_original_mu[mu][k][l] - a[mu]);
                        }
                        omega_matriz[i][j][k][l] = (1.0 / (N * N)) * sum_omega;
                    }
                }
            }
        }
    }
    return omega_matriz;
}

// Función para aplicar el algoritmo de Metropolis
void algoritmo_Metropolis(vector < vector < int > > s, vector < vector < vector < vector < double > > > > omega_matriz, int N, mt19937_64& generator, ofstream& fichero_neuronas)
{
    //PASO 0: Asigno un valor inicial a la temperatura y creo una configuración inicial de la red de neuronas

    // Asigno un valor inicial a la temperatura
    double T = 0.0001;

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