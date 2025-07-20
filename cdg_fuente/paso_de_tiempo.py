import galerkin_discontinuo
import numpy as np

def calcular_vector_residual(h__, u__, matriz_de_rigidez__):
    
    # computing stiffness vectors
    stiffness_vector_1_, stiffness_vector_2_ = galerkin_discontinuo.compute_stiffness_vectors(h__, u__, matriz_de_rigidez__)

    # computing numerical flux
    numerical_flux_vector_1_, numerical_flux_vector_2_ = galerkin_discontinuo.compute_numerical_flux_vectors(h__, u__)

    # computing residual vector
    residual_vector_1_ = stiffness_vector_1_ - numerical_flux_vector_1_
    residual_vector_2_ = stiffness_vector_2_ - numerical_flux_vector_2_

    return residual_vector_1_, residual_vector_2_

def compute_dhdt_du_dt(_h, _u, _matriz_de_rigidez, _matriz_de_masa_inversa):

    _N_elementos, _N_nodos = _h.shape

    vector_residual_1, vector_residual_2 = calcular_vector_residual(_h, _u, _matriz_de_rigidez)

    d1U_dt = np.zeros((_N_elementos, _N_nodos))
    d2U_dt = np.zeros((_N_elementos, _N_nodos))

    for i in range(_N_elementos):
        d1U_dt[i] = _matriz_de_masa_inversa @ vector_residual_1[i]
        d2U_dt[i] = _matriz_de_masa_inversa @ vector_residual_2[i]

    # inicializando los valores de dh/dt y du/dt
    dh_dt = np.zeros((_N_elementos, _N_nodos))
    du_dt = np.zeros((_N_elementos, _N_nodos))

    # calculando la derivada temporal de h y u
    # ... Escribe aqui dh_dt y du_dt ...
    # ... recuerda que 1U=h y 2U=hu por lo que du_dt = ( d2U_dt - u * dh_dt ) / h ...
    dh_dt = d1U_dt
    du_dt = np.where( _h == 0 , 0 , ( d2U_dt - _u * dh_dt ) / _h )

    return dh_dt, du_dt
