import numpy as np

def lagrange_basis_derivative(nodes, i, x):
    """
    Compute the derivative of the Lagrange basis function for a given node and value of x.

    Parameters:
        nodes (array-like): Array of nodes in the element.
        i (int): Index of the node for which to compute the derivative.
        x (float): Value of x for which to compute the derivative.

    Returns:
        float: The derivative of the Lagrange basis function at the given node and value of x.
    """
    basis_derivative = 0
    for j in range(len(nodes)):
        if j != i:
            pc = 1
            for k in range(len(nodes)):
                if k != j and k != i:
                    pc *= (x-nodes[k])/(nodes[i]-nodes[k])
            basis_derivative += pc/(nodes[i]-nodes[j])
    return basis_derivative

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
    for j in range(len(nodos)):
        # si j es diferente de i
        if j != i:
            # multiplica el polinomio de lagrange i evaluado en x por el factor correspondiente
            polinomio_de_lagrange_i_evaluado_en_x *= (x - nodos[j]) / (nodos[i] - nodos[j])

    #-----------------------------------------------------------------------------------------

    return polinomio_de_lagrange_i_evaluado_en_x

def generate_reference_space(malla_, cuadratura_de_gauss_):

    # guardando la cuadratura de Gauss en el espacio físico
    cuadratura_de_gauss_espacio_fisico = 0.5 * ( malla_[0][-1] - malla_[0][0] ) * cuadratura_de_gauss_ + 0.5 * ( malla_[0][-1] + malla_[0][0])

    # evaluando la función base en los puntos de cuadratura de Gauss
    polinomios_de_lagrange_en_cuadratura_de_gauss_ = np.array([[polinomios_de_lagrange(malla_[0], i, x) for x in cuadratura_de_gauss_espacio_fisico] for i in range(len(malla_[0]))])

    # evaluando la derivada en x de la función base evaluada en los puntos de cuadratura de Gauss
    derivada_x_polinomios_de_lagrange_en_cuadratura_de_gauss_ = np.array([[lagrange_basis_derivative(malla_[0], i, x) for x in cuadratura_de_gauss_espacio_fisico] for i in range(len(malla_[0]))])

    return polinomios_de_lagrange_en_cuadratura_de_gauss_, derivada_x_polinomios_de_lagrange_en_cuadratura_de_gauss_