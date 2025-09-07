import subprocess
import json
import sys
import random
from Bio import AlignIO, SeqIO
from statistics import mode, StatisticsError

# CARGA DE CONFIGURACIÓN

def cargar_configuracion():
    """Carga y valida la configuración desde parametros.json."""
    try:
        with open("parametros.json", "r", encoding="utf-8") as archivo:
            configuracion = json.load(archivo)

        claves_requeridas = {
            "filtro": ["archivo_salida"],
            "mafft": ["hilos", "ep", "op", "salida"],
            "ugene": ["umbral", "formato", "archivo_salida"],
            "biopython_consensus": ["umbral", "habilitado", "archivo_salida"]
        }
        
        for seccion, claves in claves_requeridas.items():
            if not all(clave in configuracion.get(seccion, {}) for clave in claves):
                faltantes = [c for c in claves if c not in configuracion[seccion]]
                raise KeyError(f"Faltan claves en 'parametros.json': {faltantes}")

        return configuracion

    except FileNotFoundError:
        print("❌ Error: Archivo 'parametros.json' no encontrado.")
        sys.exit(1)
    except json.JSONDecodeError:
        print("❌ Error: 'parametros.json' no es un JSON válido.")
        sys.exit(1)
    except KeyError as e:
        print(f"❌ Error en configuración: {str(e)}")
        sys.exit(1)

# -------------------------------------------------
# EJECUCIÓN DE COMANDOS
# -------------------------------------------------
def ejecutar_comando(comando):
    """Ejecuta un comando en la terminal y maneja errores."""
    try:
        print(f"🔹 Ejecutando: {comando}")
        subprocess.run(comando, check=True, shell=True, capture_output=True, text=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Error al ejecutar: {e}")
        print(f"   Stderr: {e.stderr}")
        return False

def ejecutar_comando_ugene(config, archivo_entrada):
    """Ejecuta UGENE con algoritmo de consenso."""
    salida_ugene = config["ugene"]["archivo_salida"]

    comando = (
        f'ugene --task=extract_consensus_sequence '
        f'--in={archivo_entrada} '
        f'--out={salida_ugene} '
        f'--format={config["ugene"]["formato"]} '
        f'--keep-gaps={str(config["ugene"].get("mantener_gaps", False)).lower()} '
        f'--threshold={config["ugene"]["umbral"]}'
    )
    return ejecutar_comando(comando)

# -------------------------------------------------
# FUNCIÓN DE RECORTE
# -------------------------------------------------
def es_codon_valido(secuencia, posicion):
    """Verifica que la posición no tenga gaps en el triplete"""
    return (posicion + 2 < len(secuencia) and 
            secuencia[posicion] != '-' and 
            secuencia[posicion+1] != '-' and 
            secuencia[posicion+2] != '-')

def recortar_secuencias(archivo_entrada, archivo_salida, config):
    """
    Recorta secuencias usando posiciones fijas (si están definidas) o codones (si no lo están).
    """
    try:
        params = config["mafft"]["procesar_codones"]
        registros = list(SeqIO.parse(archivo_entrada, "fasta"))
        
        if not registros:
            print("❌ Error: No se encontraron secuencias en el archivo de entrada")
            return False

        pos_inicio_fijo = params.get("posicion_inicio_fijo")
        pos_fin_fijo = params.get("posicion_fin_fijo")
        codon_inicio = params["codon_inicio"][0]
        codones_parada = set(params["codones_parada"])

        posiciones_inicio, posiciones_fin = [], []

        for registro in registros:
            secuencia = str(registro.seq).upper()
            if len(secuencia) < 3:
                print(f"⚠️  Secuencia {registro.id} demasiado corta")
                continue

            # --- INICIO ---
            if pos_inicio_fijo is None:
                pos_inicio = -1
                for i in range(0, len(secuencia)-2):
                    if es_codon_valido(secuencia, i) and secuencia[i:i+3] == codon_inicio:
                        pos_inicio = i
                        break
                posiciones_inicio.append(pos_inicio if pos_inicio != -1 else 0)

            # --- FIN ---
            if pos_fin_fijo is None:
                pos_fin = -1
                for i in range(len(secuencia)-3, 0, -1):
                    if es_codon_valido(secuencia, i) and secuencia[i:i+3] in codones_parada:
                        pos_fin = i + 2
                        break
                posiciones_fin.append(pos_fin if pos_fin != -1 else len(secuencia)-1)

        try:
            inicio_comun = int(pos_inicio_fijo) if pos_inicio_fijo is not None else mode(posiciones_inicio)
        except StatisticsError:
            inicio_comun = min(posiciones_inicio) if posiciones_inicio else 0
            
        try:
            fin_comun = int(pos_fin_fijo-2) if pos_fin_fijo is not None else mode(posiciones_fin)
        except StatisticsError:
            fin_comun = max(posiciones_fin) if posiciones_fin else len(secuencia)-1

        if inicio_comun < 0 or fin_comun < 0 or inicio_comun >= fin_comun:
            print("❌ Error en las posiciones de corte")
            return False

        with open(archivo_salida, "w") as salida:
            for registro in registros:
                secuencia = str(registro.seq)
                fin = min(fin_comun, len(secuencia)-1)
                recortada = secuencia[inicio_comun:fin + 1]
                salida.write(f">{registro.id}\n{recortada}\n")

        print(f"✅ Secuencias recortadas guardadas en: {archivo_salida}")
        return True

    except Exception as e:
        print(f"❌ Error inesperado al recortar: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

# -------------------------------------------------
# FUNCIÓN DE CONSENSO LEVITSKY
# -------------------------------------------------
def generar_consenso_levitsky(archivo_entrada, archivo_salida, config):
    """Genera consenso estilo Levitsky con IUPAC."""
    try:
        params = config["biopython_consensus"]
        umbral = params.get("umbral", 0.6)
        ignorar_gaps = params.get("ignorar_gaps", True)

        print(f"🔬 Generando consenso IUPAC (Levitsky, umbral={umbral})...")

        alineamiento = AlignIO.read(archivo_entrada, "fasta")
        longitud = alineamiento.get_alignment_length()
        n_secuencias = len(alineamiento)
        if n_secuencias == 0:
            print("❌ Error: Alineamiento vacío")
            return False

        iupac_map = {
            'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
            'AC': 'M', 'AG': 'R', 'AT': 'W',
            'CG': 'S', 'CT': 'Y', 'GT': 'K',
            'ACG': 'V', 'ACT': 'H', 'AGT': 'D', 'CGT': 'B',
            'ACGT': 'N'
        }

        consenso = []
        for pos in range(longitud):
            columna = alineamiento[:, pos]
            conteo = {'A':0, 'C':0, 'G':0, 'T':0, '-':0, 'N':0}
            for nt in columna:
                nt = nt.upper()
                if nt in conteo:
                    conteo[nt] += 1
                else:
                    conteo['N'] += 1

            total_valido = n_secuencias
            if ignorar_gaps:
                total_valido -= conteo['-']

            if total_valido == 0:
                consenso.append('-')
                continue

            perfil = {base: conteo[base] / total_valido for base in ['A','C','G','T']}
            bases_consenso = [b for b,f in perfil.items() if f >= umbral]

            if not bases_consenso:
                sorted_bases = sorted(perfil.items(), key=lambda x: x[1], reverse=True)
                top_bases = [b for b,f in sorted_bases if f > 0][:3]
                clave = ''.join(sorted(top_bases))
                consenso.append(iupac_map.get(clave, 'N'))
            elif len(bases_consenso) == 1:
                consenso.append(bases_consenso[0])
            else:
                clave = ''.join(sorted(bases_consenso))
                consenso.append(iupac_map.get(clave, 'N'))

        with open(archivo_salida, "w") as salida:
            salida.write(f">consenso_levitsky_umbral_{umbral}\n")
            salida.write(''.join(consenso) + "\n")

        print(f"✅ Consenso Levitsky generado ({len(consenso)} bp)")
        return True

    except Exception as e:
        print(f"❌ Error en consenso Levitsky: {str(e)}")
        return False

# -------------------------------------------------
# FLUJO DE TRABAJO (MAIN)
# -------------------------------------------------
if __name__ == "__main__":
    config = cargar_configuracion()

    alineamiento_mafft = "alineamiento_MAFFT.fa"
    alineamiento_procesado = "alineamiento_procesado.fa"

    # Paso 1: Alineamiento con MAFFT, puede reemplazar en parametros "auto" por "genafpair", "localpair" o "globalpair"
    comando_mafft = (
        f'mafft --{config["mafft"]["metodo"]} --ep {config["mafft"]["ep"]} {config["mafft"]["opcionales"]} '
        f'--op {config["mafft"]["op"]} --thread {config["mafft"]["hilos"]} '
        f'--out {alineamiento_mafft} {config["filtro"]["archivo_salida"]}'
        
    )
    if not ejecutar_comando(comando_mafft):
        sys.exit(1)

    # Paso 2: Recorte de secuencias
    if not recortar_secuencias(alineamiento_mafft, alineamiento_procesado, config):
        sys.exit(1)

    # Paso 3: Consenso con UGENE
    consenso_ugene_ok = ejecutar_comando_ugene(config, alineamiento_procesado)

    # Paso 4: Consenso con Biopython (Levitsky)
    consenso_levitsky_ok = False
    if config["biopython_consensus"]["habilitado"]:
        consenso_levitsky_ok = generar_consenso_levitsky(
            alineamiento_procesado, 
            config["biopython_consensus"]["archivo_salida"], 
            config
        )

    # 📌 Reporte final
    print("\n✅ Flujo de trabajo completado. Archivos generados:")
    print(f"- Alineamiento MAFFT (original): {alineamiento_mafft}")
    print(f"- Alineamiento procesado (recortado): {alineamiento_procesado}")

    if consenso_ugene_ok:
        print(f"✅ Consenso UGENE: {config['ugene']['archivo_salida']}")
    else:
        print("⚠️  Consenso UGENE no se generó.")

    if config["biopython_consensus"]["habilitado"]:
        if consenso_levitsky_ok:
            print(f"✅ Consenso Levitsky (Biopython): {config['biopython_consensus']['archivo_salida']}")
        else:
            print("⚠️  Consenso Levitsky no se generó.")