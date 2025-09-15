#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Batch para calcular los valores Gamma para todas las regiones en las carpetas Data/Clouds y Data/Holes.

Desarrollado por:
    Dr. Gerardo Tinoco Guerrero
    Universidad Michoacana de San Nicolás de Hidalgo
    gerardo.tinoco@umich.mx

Fecha: Mayo, 2024
"""

import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Scripts.Gammas import Gammas
import time

# Definir el operador diferencial (Laplaciano)
L = np.array([[0], [0], [1], [0], [1]])

# Obtener el directorio base del proyecto
script_dir = os.path.dirname(os.path.abspath(__file__))
base_dir = os.path.dirname(script_dir)

# Crear directorio para resultados si no existe
results_dir = os.path.join(base_dir, "Results/Gammas")
os.makedirs(results_dir, exist_ok=True)

# Función para procesar una región
def process_region(data_dir, region_name, neighbors_count=8):
    print(f"Procesando región: {region_name} en {data_dir}")
    
    # Cargar datos de coordenadas
    p_file = os.path.join(data_dir, f"{region_name}_p.csv")
    p = np.loadtxt(p_file, delimiter=',')  # Formato: x, y, flag
    
    # Cargar datos de vecinos
    vec_file = os.path.join(data_dir, f"{region_name}_neighbors_{neighbors_count}.csv")
    vec = np.loadtxt(vec_file, delimiter=',', dtype=int)  # Índices de vecinos
    
    # Calcular Gammas
    start_time = time.time()
    gamma_values = Gammas(p, vec, L)
    end_time = time.time()
    
    # Verificar si el cálculo fue exitoso
    if gamma_values is not None:
        # Guardar resultados
        output_file = os.path.join(results_dir, f"{region_name}_gammas_{neighbors_count}.csv")
        np.savetxt(output_file, gamma_values, delimiter=',')
        
        # Generar visualización
        plt.figure(figsize=(10, 8))
        plt.scatter(p[:, 0], p[:, 1], c='blue', s=1)
        plt.title(f"Región {region_name} - Gammas calculados")
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.axis('equal')
        plt.savefig(os.path.join(results_dir, f"{region_name}_gammas_{neighbors_count}.png"), dpi=300)
        plt.close()
        
        print(f"  ✓ Gammas calculados y guardados en {output_file}")
        print(f"  ✓ Tiempo de cálculo: {end_time - start_time:.2f} segundos")
        return True
    else:
        print(f"  ✗ Error en el cálculo de Gammas para {region_name}")
        return False

# Función principal
def main():
    # Definir directorios de datos usando rutas absolutas
    # Obtener el directorio base del proyecto
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.dirname(script_dir)
    
    data_dirs = [
        os.path.join(base_dir, "Data/Clouds"),
        os.path.join(base_dir, "Data/Holes")
    ]
    
    # Contador de resultados
    successful = 0
    failed = 0
    
    # Procesar cada directorio
    for data_dir in data_dirs:
        print(f"\nProcesando directorio: {data_dir}")
        
        # Obtener todas las regiones únicas en el directorio
        regions = set()
        for file in os.listdir(data_dir):
            if file.endswith("_p.csv"):
                regions.add(file.split("_p.csv")[0])
        
        # Procesar cada región
        for region in sorted(regions):
            if process_region(data_dir, region):
                successful += 1
            else:
                failed += 1
    
    # Resumen final
    print("\n" + "="*50)
    print(f"Procesamiento completado:")
    print(f"  - Regiones procesadas exitosamente: {successful}")
    print(f"  - Regiones con errores: {failed}")
    print("="*50)

# Ejecutar el programa principal
if __name__ == "__main__":
    main()