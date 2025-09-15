import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import time
import psutil
import numpy as np
import pandas as pd
import json
from datetime import datetime
import Errors as Errors
from mGFD import Stationary, TimeDerivative1, TimeDerivative2

def measure_performance(func, *args, **kwargs):
    """
    Mide tiempo de ejecución y uso de memoria de una función
    """
    # Obtener memoria inicial
    process = psutil.Process()
    memory_before = process.memory_info().rss / 1024 / 1024  # MB
    
    # Medir tiempo de ejecución
    start_time = time.time()
    result = func(*args, **kwargs)
    end_time = time.time()
    
    # Obtener memoria final
    memory_after = process.memory_info().rss / 1024 / 1024  # MB
    
    execution_time = end_time - start_time
    memory_used = memory_after - memory_before
    peak_memory = memory_after
    
    return result, execution_time, memory_used, peak_memory

def benchmark_poisson_equation(region, data_type):
    """
    Benchmark para la ecuación de Poisson
    """
    # Cargar datos
    data_path = f'Data/{data_type}/'
    p_file = os.path.join(data_path, f'{region}_p.csv')
    tt_file = os.path.join(data_path, f'{region}_tt.csv')
    
    if not (os.path.exists(p_file) and os.path.exists(tt_file)):
        return None
    
    p = np.genfromtxt(p_file, delimiter=',', skip_header=0)
    tt = np.genfromtxt(tt_file, delimiter=',', skip_header=0)
    
    # Parámetros del problema (exactamente como en run_Poisson.py)
    phi = lambda x, y: 2*np.exp(2*x + y)  # Condición de frontera
    f = lambda x, y: 10*np.exp(2*x + y)   # Lado derecho de la ecuación
    L = np.vstack([[0], [0], [2], [0], [2], [0]])  # Operador
    
    # Resolver con medición de rendimiento
    (u_ap, u_ex, vec), exec_time, memory_used, peak_memory = measure_performance(
        Stationary, p, phi, f, operator=L, triangulation=False, tt=None
    )
    
    # Calcular error
    er = Errors.Cloud_Stationary(p, vec, u_ap, u_ex)
    avg_error = np.mean(er)
    
    return {
        'equation': 'Poisson',
        'region': region,
        'data_type': data_type,
        'num_points': len(p),
        'boundary_condition': 'phi = 2*exp(2*x + y)',
        'rhs_function': 'f = 10*exp(2*x + y)',
        'operator': 'L = [0, 0, 2, 0, 2, 0]',
        'execution_time_seconds': exec_time,
        'memory_used_mb': memory_used,
        'peak_memory_mb': peak_memory,
        'avg_numerical_error': avg_error
    }

def benchmark_heat_equation(region, data_type):
    """
    Benchmark para la ecuación de calor
    """
    # Cargar datos
    data_path = f'Data/{data_type}/'
    p_file = os.path.join(data_path, f'{region}_p.csv')
    tt_file = os.path.join(data_path, f'{region}_tt.csv')
    
    if not (os.path.exists(p_file) and os.path.exists(tt_file)):
        return None
    
    p = np.genfromtxt(p_file, delimiter=',', skip_header=0)
    tt = np.genfromtxt(tt_file, delimiter=',', skip_header=0)
    
    # Parámetros del problema (exactamente como en run_Heat.py)
    v = 0.2  # Coeficiente de difusión
    t = 2000  # Número de pasos de tiempo
    f = lambda x, y, t, coef: np.exp(-2*np.pi**2*coef[0]*t)*np.cos(np.pi*x)*np.cos(np.pi*y)
    L = np.vstack([[0], [0], [2*v], [0], [2*v], [0]])  # Operador
    
    # Resolver con medición de rendimiento
    (u_ap, u_ex, vec), exec_time, memory_used, peak_memory = measure_performance(
        TimeDerivative1, p, f, t, [v], operator=L, triangulation=False, tt=[], implicit=False, lam=0.5
    )
    
    # Calcular error
    er = Errors.Cloud_Transient(p, vec, u_ap, u_ex)
    avg_error = np.mean(er)
    
    return {
        'equation': 'Heat',
        'region': region,
        'data_type': data_type,
        'num_points': len(p),
        'diffusion_coefficient': v,
        'time_steps': t,
        'function': 'f = exp(-2*pi^2*v*t)*cos(pi*x)*cos(pi*y)',
        'operator': f'L = [0, 0, {2*v}, 0, {2*v}, 0]',
        'implicit': False,
        'lambda': 0.5,
        'execution_time_seconds': exec_time,
        'memory_used_mb': memory_used,
        'peak_memory_mb': peak_memory,
        'time_per_step_ms': (exec_time * 1000) / t,
        'avg_numerical_error': avg_error
    }

def benchmark_advdif_equation(region, data_type):
    """
    Benchmark para la ecuación de Advección-Difusión
    """
    # Cargar datos
    data_path = f'Data/{data_type}/'
    p_file = os.path.join(data_path, f'{region}_p.csv')
    tt_file = os.path.join(data_path, f'{region}_tt.csv')
    
    if not (os.path.exists(p_file) and os.path.exists(tt_file)):
        return None
    
    p = np.genfromtxt(p_file, delimiter=',', skip_header=0)
    tt = np.genfromtxt(tt_file, delimiter=',', skip_header=0)
    
    # Parámetros del problema (exactamente como en run_AdvDif.py)
    v = 0.1  # Coeficiente de difusión
    a = 0.3  # Velocidad de transporte en dirección x
    b = 0.2  # Velocidad de transporte en dirección y
    t = 2000  # Número de pasos de tiempo
    f = lambda x, y, t, coef: (1/(4*t+1))*np.exp(-(x-coef[1]*t-0.5)**2/(coef[0]*(4*t+1)) - (y-coef[2]*t-0.5)**2/(coef[0]*(4*t+1)))
    L = np.vstack([[-a], [-b], [2*v], [0], [2*v], [0]])  # Operador
    
    # Resolver con medición de rendimiento
    (u_ap, u_ex, vec), exec_time, memory_used, peak_memory = measure_performance(
        TimeDerivative1, p, f, t, [v, a, b], operator=L, triangulation=False, tt=[], implicit=False, lam=0.5
    )
    
    # Calcular error
    er = Errors.Cloud_Transient(p, vec, u_ap, u_ex)
    avg_error = np.mean(er)
    
    return {
        'equation': 'Advection-Diffusion',
        'region': region,
        'data_type': data_type,
        'num_points': len(p),
        'diffusion_coefficient': v,
        'transport_velocity_x': a,
        'transport_velocity_y': b,
        'time_steps': t,
        'function': 'f = (1/(4*t+1))*exp(-(x-a*t-0.5)^2/(v*(4*t+1)) - (y-b*t-0.5)^2/(v*(4*t+1)))',
        'operator': f'L = [{-a}, {-b}, {2*v}, 0, {2*v}, 0]',
        'implicit': False,
        'lambda': 0.5,
        'execution_time_seconds': exec_time,
        'memory_used_mb': memory_used,
        'peak_memory_mb': peak_memory,
        'time_per_step_ms': (exec_time * 1000) / t,
        'avg_numerical_error': avg_error
    }

def benchmark_wave_equation(region, data_type):
    """
    Benchmark para la ecuación de onda
    """
    # Cargar datos
    data_path = f'Data/{data_type}/'
    p_file = os.path.join(data_path, f'{region}_p.csv')
    tt_file = os.path.join(data_path, f'{region}_tt.csv')
    
    if not (os.path.exists(p_file) and os.path.exists(tt_file)):
        return None
    
    p = np.genfromtxt(p_file, delimiter=',', skip_header=0)
    tt = np.genfromtxt(tt_file, delimiter=',', skip_header=0)
    
    # Parámetros del problema (exactamente como en run_Wave.py)
    c = np.sqrt(1/2)  # Coeficiente de propagación de onda
    t = 2000  # Número de pasos de tiempo
    f = lambda x, y, t, coef: np.cos(np.pi*t)*np.sin(np.pi*(x+y))
    g = lambda x, y, t, coef: -np.pi*np.sin(np.pi*t)*np.sin(np.pi*(x+y))
    L = np.vstack([[0], [0], [2*c**2], [0], [2*c**2], [0]])  # Operador
    
    # Resolver con medición de rendimiento
    (u_ap, u_ex, vec), exec_time, memory_used, peak_memory = measure_performance(
        TimeDerivative2, p, f, g, t, [c], operator=L, triangulation=False, tt=None, implicit=False, lam=1
    )
    
    # Calcular error
    er = Errors.Cloud_Transient(p, vec, u_ap, u_ex)
    avg_error = np.mean(er)
    
    return {
        'equation': 'Wave',
        'region': region,
        'data_type': data_type,
        'num_points': len(p),
        'wave_coefficient': c,
        'time_steps': t,
        'function_f': 'f = cos(pi*t)*sin(pi*(x+y))',
        'function_g': 'g = -pi*sin(pi*t)*sin(pi*(x+y))',
        'operator': f'L = [0, 0, {2*c**2}, 0, {2*c**2}, 0]',
        'implicit': False,
        'lambda': 1,
        'execution_time_seconds': exec_time,
        'memory_used_mb': memory_used,
        'peak_memory_mb': peak_memory,
        'time_per_step_ms': (exec_time * 1000) / t,
        'avg_numerical_error': avg_error
    }

def run_comprehensive_benchmark(regions=None, data_types=None, equations=None, save_results=True):
    """
    Ejecuta el benchmark completo
    """
    if regions is None:
        regions = ['BAN', 'BLU', 'CUA', 'CUI', 'ENG', 'GIB', 'HAB', 'MIC', 'PAT', 'TIT', 'TOB', 'UCH', 'VAL', 'ZIR']
    
    if data_types is None:
        data_types = ['Clouds', 'Holes']
    
    if equations is None:
        equations = ['Poisson', 'Heat', 'AdvDif', 'Wave']
    
    # Mapeo de ecuaciones a funciones
    equation_functions = {
        'Poisson': benchmark_poisson_equation,
        'Heat': benchmark_heat_equation,
        'AdvDif': benchmark_advdif_equation,
        'Wave': benchmark_wave_equation
    }
    
    results = []
    
    print("Iniciando benchmark completo...")
    print(f"Regiones: {regions}")
    print(f"Tipos de datos: {data_types}")
    print(f"Ecuaciones: {equations}")
    
    total_combinations = len(regions) * len(data_types) * len(equations)
    current = 0
    
    # Tiempo total del benchmark
    benchmark_start_time = time.time()
    
    for data_type in data_types:
        for region in regions:
            for equation in equations:
                current += 1
                print(f"\nProcesando {current}/{total_combinations}: {region} ({data_type}) - {equation}")
                
                try:
                    benchmark_func = equation_functions[equation]
                    result = benchmark_func(region, data_type)
                    if result:
                        results.append(result)
                        print(f"  Error promedio: {result['avg_numerical_error']:.2e}")
                        print(f"  Tiempo: {result['execution_time_seconds']:.3f}s")
                        print(f"  Memoria: {result['memory_used_mb']:.1f}MB (pico: {result['peak_memory_mb']:.1f}MB)")
                        if 'time_per_step_ms' in result:
                            print(f"  Tiempo por paso: {result['time_per_step_ms']:.3f}ms")
                    else:
                        print(f"  Archivos no encontrados para {region} en {data_type}")
                except Exception as e:
                    print(f"  Error en {equation}: {e}")
    
    benchmark_end_time = time.time()
    total_benchmark_time = benchmark_end_time - benchmark_start_time
    
    if save_results and results:
        # Crear directorio para resultados
        os.makedirs('benchmark_results', exist_ok=True)
        
        # Guardar los nombres de las ecuaciones
        equation_names = {
            'Poisson': 'Poisson',
            'Heat': 'Heat',
            'AdvDif': 'Advection-Diffusion',
            'Wave': 'Wave'
        }

        # Agregar estadísticas del benchmark
        benchmark_stats = {
            'total_benchmark_time_seconds': total_benchmark_time,
            'total_tests': len(results),
            'successful_tests': len([r for r in results if r is not None]),
            'timestamp': datetime.now().isoformat(),
            'regions_tested': regions,
            'data_types_tested': data_types,
            'equations_tested': equations
        }
        
        # Guardar como JSON
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        json_file = f'benchmark_results/benchmark_{timestamp}.json'
        
        output_data = {
            'benchmark_stats': benchmark_stats,
            'results': results
        }
        
        with open(json_file, 'w') as f:
            json.dump(output_data, f, indent=2)
        
        # Guardar como CSV
        csv_file = f'benchmark_results/benchmark_{timestamp}.csv'
        df = pd.DataFrame(results)
        df.to_csv(csv_file, index=False)
        
        # Crear resumen estadístico
        summary_file = f'benchmark_results/benchmark_summary_{timestamp}.txt'
        with open(summary_file, 'w') as f:
            f.write(f"RESUMEN DEL BENCHMARK\n")
            f.write(f"=====================\n\n")
            f.write(f"Tiempo total del benchmark: {total_benchmark_time:.2f} segundos\n")
            f.write(f"Pruebas exitosas: {len(results)} de {total_combinations}\n\n")
            
            if results:
                df = pd.DataFrame(results)
                f.write(f"ESTADÍSTICAS POR ECUACIÓN:\n")
                for eq in equations:
                    eq_name = equation_names.get(eq, eq)
                    eq_results = df[df['equation'] == eq_name]
                    if not eq_results.empty:
                        f.write(f"\n{eq_name}:\n")
                        f.write(f"  Tiempo promedio: {eq_results['execution_time_seconds'].mean():.3f}s\n")
                        f.write(f"  Memoria promedio: {eq_results['memory_used_mb'].mean():.1f}MB\n")
                        f.write(f"  Error promedio: {eq_results['avg_numerical_error'].mean():.2e}\n")
        
        print(f"\nResultados guardados en:")
        print(f"  JSON: {json_file}")
        print(f"  CSV: {csv_file}")
        print(f"  Resumen: {summary_file}")
    
    print(f"\nBenchmark completado en {total_benchmark_time:.2f} segundos")
    print(f"Total de resultados: {len(results)}")
    return results

if __name__ == "__main__":
    # Verificar que psutil esté disponible
    try:
        import psutil
    except ImportError:
        print("Error: psutil no está instalado. Instálalo con: pip install psutil")
        exit(1)
    
    # Ejecutar benchmark con todas las ecuaciones
    results = run_comprehensive_benchmark()