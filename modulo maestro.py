# run_pipeline.py
import subprocess
import sys
from pathlib import Path

def ejecutar_modulo(script, descripcion):
    """Ejecuta un módulo y maneja errores."""
    print(f"\n{'='*50}\n🔹 Ejecutando: {descripcion}\n{'='*50}")
    try:
        result = subprocess.run([sys.executable, script], check=True)
        if result.returncode == 0:
            print(f"✅ {descripcion} completado exitosamente.")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Error en {script}: {e}")
        sys.exit(1)

if __name__ == "__main__":
    # Orden de ejecución
    modulos = [
        ("1-Filtracion.py", "Filtrado de secuencias"),
        ("2-Alineamiento.py", "Alineamiento con MAFFT y Consenso con UGENE"),
        ("3-Reporte.py", "Reporte de cebadores")
    ]

    # Verificar que los archivos estan
    for script, _ in modulos:
        if not Path(script).exists():
            print(f"❌ Error: No se encontró {script}")
            sys.exit(1)

    # Ejecutar pipeline
    for script, descripcion in modulos:
        if not ejecutar_modulo(script, descripcion):
            sys.exit(1)

    print("\n" + "="*50)
    print("Modulos ejecutados exitosamente. Verifique los resultados.")
    print("="*50)