<!-- README.md -->
# Técnicas Computacionales de Fluidos en Ingeniería y Astrofísica  
**Taller InspiraSTEM 2025**

| | |
|---|---|
| **Lugar** | *[Completar]* |
| **Fechas** | *[Completar]* |
| **Hora** | *[Completar]* |
| **Instructor** | **Erick Urquilla** — eurquill@vols.utk.edu |

---

## Sobre el instructor  
Soy un estudiante salvadoreño de doctorado en Física (Bains Fellow) originario de San Juan Opico, La Libertad. Bajo la asesoría del Dr. Sherwood Richers en la Universidad de Tennessee, Knoxville, estudio fenómenos cuánticos de neutrinos en supernovas y colisiones de estrellas de neutrones. He trabajado con **CFD avanzado** en el Departamento de Ingeniería Aeroespacial de UTK y con métodos de transporte de neutrinos en astrofísica.

---

## Descripción del curso  
Este curso intensivo explora **técnicas avanzadas de modelado computacional de fluidos**. A lo largo de tres días —con un enfoque 100 % práctico— los participantes completarán componentes clave de un *solver* existente para un **fluido incompresible 1-D**:

* **Discretización espacial**: Elementos Finitos Discontinuos (*Discontinuous Galerkin*).  
* **Discretización temporal**: Diferencias Finitas (Euler hacia adelante o Runge–Kutta 4).  

Además de la aplicación directa en ingeniería, los métodos se conectarán con procesos astrofísicos como el transporte de materia y fotones.

---

## Perfil de estudiantes  
| Requisito | Nivel recomendado |
|-----------|-------------------|
| Álgebra, geometría, álgebra lineal | Básico–intermedio |
| Cálculo & EDP (quinto semestre) | Intermedio |
| Conceptos básicos de fluidos | Deseable |
| Python (variables, funciones, *NumPy*, *Matplotlib*) | Sólido |

---

## Objetivos generales  
1. Comprender la formulación numérica de la ecuación de agua poco profunda 1-D.  
2. Implementar un solver DG + RK4/Euler en **Python**.  
3. Emplear integración de Gauss y bases de **Lagrange**.  
4. Analizar estabilidad y estimar errores numéricos.  
5. Transferir las técnicas a problemas astrofísicos de mayor complejidad.

---

## Plan diario y syllabus detallado  

| Día | Objetivos clave | Contenidos & actividades | Software |
|-----|-----------------|--------------------------|----------|
| **1** | • Rol del CFD en ciencia e ingeniería<br>• Generación de cuadrículas<br>• Condiciones iniciales | 1. Presentación del curso <br>2. Proyecto de grupo: *“Agua en un guacal 1-D”*<br>3. **Ejercicio 1** — función Python para generar cuadrículas<br>4. **Ejercicio 2** — codificar condiciones iniciales | Jupyter, NumPy, Matplotlib, IPython, `time` |
| **2** | • Bases de Lagrange y aproximación funcional<br>• Forma débil 1-D<br>• Integración numérica (Gauss) | 1. Galerkin continuo vs discontinuo<br>2. **Ejercicio 3** — polinomio de Lagrange evaluado en `x`<br>3. Forma débil de las ecuaciones<br>4. **Ejercicio 4** — matriz de masa por cuadratura de Gauss | Jupyter, NumPy, Matplotlib |
| **3** | • Solución temporal con RK4 / Newton<br>• Forma débil completa<br>• Cierre del proyecto | 1. Forma débil (parte II)<br>2. Método de Newton; RK4<br>3. **Ejercicio 5** — implementar Euler o RK4 para la evolución temporal | Jupyter, NumPy, Matplotlib |

> **Nota:** Cada día combina micro-exposiciones (~20 min) + bloques de codificación supervisada (~90 min) + debrief (~30 min).

---

## Estructura del repositorio  

| Carpeta / archivo | Descripción |
|-------------------|-------------|
| `notebooks/` | Cuadernos paso a paso (Día 1-3). |
| `src/` | Módulos Python del solver DG. |
| `data/` | Resultados y conjuntos de prueba. |
| `docs/` | Diapositivas y material de referencia. |
| `requirements.txt` | Dependencias mínimas. |

---

## Instalación rápida

```bash
git clone https://github.com/InspiraSTEM/taller-cfd-2025.git
cd taller-cfd-2025

python -m venv .venv
source .venv/bin/activate        # Linux/Mac
# .venv\Scripts\activate.bat     # Windows

pip install -r requirements.txt
jupyter lab notebooks/
