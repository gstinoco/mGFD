import os
import sys

# Agregamos la ruta base para poder importar los módulos
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(BASE_DIR)
sys.path.append(os.path.dirname(BASE_DIR)) # Para que mGFD sea accesible

from mGFD.io.io import load_points
import benchmark as bmk

def test_single_cloud():
    """
    Prueba rápida de todos los solvers sobre una única nube de puntos.
    """
    cloud_path = os.path.join(os.path.dirname(BASE_DIR), 'Data', 'Titicaca', '1', 'Titicaca_cloud.csv')
    if not os.path.exists(cloud_path):
        print(f"No se encontró la nube de prueba en: {cloud_path}")
        return

    print(f"Iniciando prueba rápida en una única nube: {cloud_path}")
    p = load_points(cloud_path)
    print(f"Nube cargada con {len(p)} puntos.")
    
    # Lista de funciones a probar
    tests = [
        ('Poisson', bmk.benchmark_poisson_equation),
        ('Heat (t=2000)', bmk.benchmark_heat_equation),
        ('Advection', bmk.benchmark_advection_equation),
        ('Advection-Diffusion', bmk.benchmark_advdif_equation),
        ('Wave', bmk.benchmark_wave_equation)
    ]
    
    for name, func in tests:
        print(f"\n[{name}] Ejecutando...")
        try:
            res = func(p)
            print(f"  [EXITO] Tiempo: {res['execution_time_seconds']:.3f}s | RMSE Promedio: {res['avg_numerical_error']:.3e}")
        except Exception as e:
            print(f"  [ERROR] Falló la prueba {name}: {e}")

if __name__ == '__main__':
    test_single_cloud()
