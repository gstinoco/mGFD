import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Scripts.Neighbors import Cloud
from Scripts.Neighbors import Triangulation


def generate_neighbor_files(nvec=10):
    """Genera archivos de vecinos para todas las regiones"""
    regions = ['BAN', 'BLU', 'CUA', 'CUI', 'ENG', 'GIB', 'HAB', 'MIC', 'PAT', 'TIT', 'TOB', 'UCH', 'VAL', 'ZIR']
    
    for region in regions:
        for data_type in ['Clouds', 'Holes']:
            p_file            = f'Data/{data_type}/{region}_p.csv'
            tt_file           = f'Data/{data_type}/{region}_tt.csv'
            neighbors_file    = f'Data/{data_type}/{region}_neighbors_{nvec}.csv'
            neighbors_file_tt = f'Data/{data_type}/{region}_neighbors_tt_{nvec}.csv'
            
            if os.path.exists(p_file):
                print(f"Procesando {region} ({data_type})...")
                p_data = np.loadtxt(p_file, delimiter=',')
                tt     = np.loadtxt(tt_file, delimiter=',')
                
                # Calcular vecinos usando la función Cloud
                vec_data    = Cloud(p_data, nvec)
                vec_data_tt = Triangulation(p_data, tt, nvec)
                
                # Guardar vecinos
                np.savetxt(neighbors_file, vec_data, delimiter=',', fmt='%d')
                print(f"  Vecinos guardados: {neighbors_file}")
                np.savetxt(neighbors_file_tt, vec_data_tt, delimiter=',', fmt='%d')
                print(f"  Vecinos guardados: {neighbors_file_tt}")

if __name__ == "__main__":
    generate_neighbor_files(nvec=10)