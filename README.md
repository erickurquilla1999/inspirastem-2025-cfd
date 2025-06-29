<!-- README.md -->
# Técnicas Computacionales de Fluidos en Ingeniería y Astrofísica  

**Erick Urquilla**\
PhD Student\
Department of Physics and Astronomy\
University of Tennessee, Knoxville\
eurquill@vols.utk.edu\
[sites.google.com/view/erickurquilla](https://sites.google.com/view/erickurquilla)

Soy originario de San Juan Opico, departamento de La Libertad, El Salvador. Obtuve la licenciatura en Física en la Universidad de El Salvador, donde inicié mi trayectoria investigadora junto al Dr. Raúl Ortiz, realizando simulaciones cosmológicas del Universo que incluían agujeros negros primordiales como modelos de materia oscura. Actualmente curso el doctorado en Física como Bains Fellow en la Universidad de Tennessee, Knoxville, bajo la dirección del Dr. Sherwood Richers, centrando mi trabajo en los fenómenos cuánticos de los neutrinos en supernovas y colisiones de estrellas de neutrones. Además, en el Departamento de Ingeniería Aeroespacial de UTK he elaborado proyectos de dinámica de fluidos computacional avanzada y he desarrollado métodos de transporte de neutrinos aplicados a la astrofísica.

---

## Descripción del curso  

Este curso aborda técnicas avanzadas de modelado computacional de fluidos. Con un enfoque práctico, los estudiantes trabajarán en simulaciones de un fluido incompresible en una dimensión, completando partes clave de un código existente. Se utilizarán elementos finitos (Discontinuous Galerkin) para la discretización espacial y diferencias finitas (Runge-Kutta 4 o Euler hacia adelante) para la temporal. Se cubrirán temas como la generación de mallas, la forma débil de ecuaciones diferenciales, funciones base de Lagrange, solucionadores de Riemann y estimación de errores. Las técnicas aprendidas podrán aplicarse a procesos astrofísicos más complejos, como el transporte de materia y fotones.

---

## Objetivos generales  
1. Comprender la formulación numérica de la ecuación de agua poco profunda 1-D.  
2. Implementar un solver DG + RK4/Euler en **Python**.  
3. Emplear integración de Gauss y bases de **Lagrange**.  
4. Analizar estabilidad y estimar errores numéricos.  
5. Transferir las técnicas a problemas astrofísicos de mayor complejidad.

---

## Perfil de estudiantes  

Es deseable para el curso tener fundamentos básicos de álgebra, geometría, álgebra lineal, cálculo y ecuaciones diferenciales parciales a un nivel de quinto semestre de carreras de ciencias e ingenierías. Conceptos básicos de fluidos. Habilidades de programación en Python: tipos de variables, operaciones logicas, loops, funciones, paquetes. Conocimiento de arrays numpy y gráficos de matplotlib harían la resolución de ejercicios más simple.

---

## Material de estudio pre-conferencia  

| Recurso | Enlace |
|---------| --------|
| Curso de introduccion a Google Colab y Python (tipos de variables, operaciones logicas, loops, funciones, paquetes) |<https://www.youtube.com/playlist?list=PLKd7y--oK26dWXHV0Aoi5eOgoikk3Fp-J>|
| Curso de introduccion numpy arrays |<https://www.youtube.com/watch?v=fczpaVSwGCA>|
| Ecuaciones de agua profunda (Wikipedia) |<https://en.wikipedia.org/wiki/Shallow_water_equations>|

---
