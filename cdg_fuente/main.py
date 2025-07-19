import numpy as np
import matplotlib.pyplot as pl











import bases
import galerkin_discontinuo 
import paso_de_tiempo
import graficos













# Parametros de entrada

# Dominio espacial de la simulacion 
x_inicial = 0 # (m) cooordenda inicial del dominio
x_final = 10 # (m) coordenada final del dominio

# Parametros del metodo de elementos finitos
N_elementos = 6 # numero de elementos
N_nodos = 4 # numero de nodos por elemento (por simplicidad solo consideramos N_nodos >= 2)

# Dominio temporal 
n_pasos = 5 # numero de pasos temporales
t_total = 1 # (s) tiempo final

# Integracion numerica con cuadratura de Gauss
n_nodos_cuadratura_gauss = 20

# Variables utiles
longitud_elemento = (x_final - x_inicial) / N_elementos














print(f'Generando malla \nDominio físico: [{x_inicial},{x_final}] metros\nNúmero de elementos: {N_elementos}\nNodos por elemento: {N_nodos}\n')

# Generar coordenadas de los elementos y nodos en el espacio físico
malla = np.zeros((N_elementos, N_nodos)) # nunpy array que almacenara las coordenadas de los nodos de cada elemento

#-----------------------------------------------------------------------------------------
# Escribe tu solucion al ejercicio 1 a continuacion ...
malla = np.array([np.linspace(x_inicial + i * (x_final - x_inicial) / N_elementos, x_inicial + (i + 1) * (x_final - x_inicial) / N_elementos, N_nodos) for i in range(N_elementos)])      
#-----------------------------------------------------------------------------------------















# Generando condiciones iniciales
h = np.zeros((N_elementos, N_nodos)) # Height (cm)
u = np.zeros((N_elementos, N_nodos)) # Velocity (cm/s)

#-----------------------------------------------------------------------------------------
# Escribe tu solucion al ejercicio 2 a continuacion ...
h = 1.0 + 0.1 * np.exp( - ( malla - 5.0 )**2 ) # Height (cm)
u = malla * 0.0 # Velocity (cm/s)
#-----------------------------------------------------------------------------------------














def polinomios_de_lagrange(nodos, i, x):
    '''
    Dado un conjunto de nodos en un elemento esta funcion
    evalua el polinomio de Lagrange i (es decir el polinomio de 
    lagrange que es uno cuando es evaluado en la posicion del nodo i) 
    en x la posicion x

    Parámetros:
        nodos (numpy.ndarray): Array de nodos de un elemento de la malla.
            ejemplo: nodos = np.array([x0, x1, ..., xN]) donde N es el numero de nodos por elemento
        i (int): Índice del nodo en el cual el polinomio de lagrange es uno.
            ejemplo: i = 1
        x (float): Punto en el cual se evalúa la función base.
            ejemplo: x = 0.5

    Retorna:
        float: Valor de la función base de Lagrange i en el punto x.
    '''
    
    # inicializa el polinomio de lagrange i evaluado en x
    polinomio_de_lagrange_i_evaluado_en_x = 1.0

    #-----------------------------------------------------------------------------------------
    # Escribe tu solucion al ejercicio 3 a continuacion ...
    
    # calcula el polinomio de lagrange i evaluado en x
    for j in range(N_nodos):
        # si j es diferente de i
        if j != i:
            # multiplica el polinomio de lagrange i evaluado en x por el factor correspondiente
            polinomio_de_lagrange_i_evaluado_en_x *= (x - nodos[j]) / (nodos[i] - nodos[j])

    #-----------------------------------------------------------------------------------------

    return polinomio_de_lagrange_i_evaluado_en_x

# Verifica que tu solucion satisface la propiedad de los polinomios de lagrange 
# \phi_i(x_j) = \delta_{ij}
nodos_prueba = malla[0]
i_prueba = 1
x_prueba = 0
phi_prueba = polinomios_de_lagrange(nodos_prueba, i_prueba, x_prueba)
print(f"Para los nodos {nodos_prueba}, el polinomio de lagrange i={i_prueba} evaluado en x={x_prueba} es {phi_prueba}")









cuadratura_de_gauss, pesos_de_gauss = np.polynomial.legendre.leggauss(n_nodos_cuadratura_gauss)









polinomios_de_lagrange_en_cuadratura_de_gauss, derivada_x_polinomios_de_lagrange_en_cuadratura_de_gauss = bases.generate_reference_space(malla, cuadratura_de_gauss, polinomios_de_lagrange)












def calcula_inversa_matriz_de_masa(longitud_elemento_, pesos_de_gauss_, polinomios_de_lagrange_en_cuadratura_de_gauss_):
    """
    Calcula la inversa de la matriz de masa para un elemento finito unidimensional.

    Parametros:

    longitud_elemento_ (float): La longitud del elemento finito.
        ejemplo: longitud_elemento = 1.0

    pesos_de_gauss_ (numpy.ndarray): Los pesos de Gauss utilizados en la cuadratura de Gauss.
        ejemplo: pesos_de_gauss = [w0, w1, ..., wN] donde N es el numero de nodos de cuadratura de Gauss
    
    polinomios_de_lagrange_en_cuadratura_de_gauss_ (numpy.ndarray): Los valores de los polinomios de Lagrange evaluados en los puntos de cuadratura de Gauss.
        ejemplo: 
            polinomios_de_lagrange_en_cuadratura_de_gauss_ = [[phi_0(x_0), phi_0(x_1), ..., phi_0(x_N)], [phi_1(x_0), phi_1(x_1), ..., phi_1(x_N)], ..., [phi_N(x_0), phi_N(x_1), ..., phi_N(x_N)]]
            polinomios_de_lagrange_en_cuadratura_de_gauss_[i] es un array con los valores de phi_i(x) evaluados en los nodos de cuadratura de Gauss.
            polinomios_de_lagrange_en_cuadratura_de_gauss_[i] = [phi_i(x_0), phi_i(x_1), ..., phi_i(x_N)]

    Retorna:
    
    numpy.ndarray: La inversa de la matriz de masa.
    """

    # Inicializa la matriz de masa
    matriz_de_masa = np.zeros((N_nodos, N_nodos))

    #-----------------------------------------------------------------------------------------
    # Escribe tu solucion al ejercicio 4 a continuacion ...

    # Calcula la matriz de masa
    for i in range(N_nodos):
        for j in range(N_nodos):
            # Calcula la integral de phi_i(x) * phi_j(x) dx
            matriz_de_masa[i, j] = 0.5 * longitud_elemento_ * np.sum(pesos_de_gauss_ * polinomios_de_lagrange_en_cuadratura_de_gauss_[i] * polinomios_de_lagrange_en_cuadratura_de_gauss_[j])
    # Calcula la inversa de la matriz de masa
    inversa_matriz_de_masa = np.linalg.inv(matriz_de_masa)

    #-----------------------------------------------------------------------------------------

    return inversa_matriz_de_masa

# calcula la matriz de masa inversa, esto es la inversa de la matriz: M_ij = integral phi_i(x) phi_j(x) dx
matriz_de_masa_inversa = calcula_inversa_matriz_de_masa(longitud_elemento, pesos_de_gauss, polinomios_de_lagrange_en_cuadratura_de_gauss)




















# calcula la matriz de rigidez S_ij = integral (dphi_i(x)/dx)phi_j(x) dx
matriz_de_rigidez = galerkin_discontinuo.calcula_matrix_de_rigidez(longitud_elemento, pesos_de_gauss, polinomios_de_lagrange_en_cuadratura_de_gauss, derivada_x_polinomios_de_lagrange_en_cuadratura_de_gauss)



















# time step
time_step = np.array(t_total/n_pasos) 

# evolving in time the PDE
for number_of_t_step in np.arange(n_pasos):

    # calculando el vector residual
    vector_residual_1, vector_residual_2 = paso_de_tiempo.calcular_vector_residual(h, u, matriz_de_rigidez)

    # -----------------------------------------------------------------------------------------
    # Escribe tu solucion al ejercicio 5 a continuacion ...

    # inicializando los valores de d1U/dt y d2U/dt
    d1U_dt = np.zeros((N_elementos, N_nodos))
    d2U_dt = np.zeros((N_elementos, N_nodos))

    # calculando la derivada temporal de 1U y 2U ...
    # ... Escribe aqui d1U_dt y d2U_dt ... 
    # ... Haz un loop sobre cada elemento ...
    # ... Utiliza la matriz de masa inversa y el vector residual en cada elemento para calcular las derivadas ...
    for i in range(N_elementos):
        d1U_dt[i] = matriz_de_masa_inversa @ vector_residual_1[i]
        d2U_dt[i] = matriz_de_masa_inversa @ vector_residual_2[i]

    # inicializando los valores de dh/dt y du/dt
    dh_dt = np.zeros((N_elementos, N_nodos))
    du_dt = np.zeros((N_elementos, N_nodos))

    # calculando la derivada temporal de h y u
    # ... Escribe aqui dh_dt y du_dt ...
    # ... recuerda que 1U=h y 2U=hu por lo que du_dt = ( d2U_dt - u * dh_dt ) / h ...
    dh_dt = d1U_dt
    du_dt = np.where( h == 0 , 0 , ( d2U_dt - u * dh_dt ) / h )

    # inicializando los valores de h y u en el siguiente paso de tiempo
    h_nuevo = np.zeros((N_elementos, N_nodos))
    u_nuevo = np.zeros((N_elementos, N_nodos))

    # calculando los valores de h y u en el siguiente paso de tiempo con el metodo de Euler explicito
    # ... Escribe aqui h_nuevo y u_nuevo dados por el metodo de Euler explicito ...
    h_nuevo = h + dh_dt * time_step
    u_nuevo = u + du_dt * time_step
    
    # -----------------------------------------------------------------------------------------

    # graficando la simulacion
    graficos.plot_simulation(malla, h_nuevo, u_nuevo, N_elementos, time_step, number_of_t_step)

    # actualizando los valores de h y u
    h = h_nuevo
    u = u_nuevo

print(f'Done')



