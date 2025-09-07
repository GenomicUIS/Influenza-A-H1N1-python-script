import re
import json

# Cargar configuración desde 'parametros.json'
with open("parametros.json", "r") as config_file:
    config = json.load(config_file)

# Parámetros desde JSON
ARCHIVO_ENTRADA = config["filtro"]["archivo_entrada"]
ARCHIVO_SALIDA = config["filtro"]["archivo_salida"]
ARCHIVO_SALIDA_N = config["filtro"].get("archivo_salida_N", "secuencias_con_N.fasta")
AÑO = config["filtro"]["periodo"]
MES = config["filtro"].get("mes")  # Puede ser null → None

def extraer_fecha(cabecera: str) -> str | None:
    """
    Extrae la fecha desde un campo delimitado por '|' en la cabecera FASTA.
    Acepta YYYY, YYYY-MM, YYYY-MM-DD con separador '-' o '/'.
    Devuelve 'YYYY-MM-DD' normalizado (rellena con 00 lo que falte).
    """
    campos = [c.strip() for c in cabecera.split('|')]
    patron = re.compile(r'^\d{4}(?:[-/]\d{1,2}){0,2}$')  # regex para fecha

    for campo in campos:
        if patron.fullmatch(campo):
            partes = re.split(r'[-/]', campo)
            año = partes[0]
            mes = partes[1].zfill(2) if len(partes) > 1 else '00'
            dia = partes[2].zfill(2) if len(partes) > 2 else '00'
            return f"{año}-{mes}-{dia}"
    return None

def validar_fecha(fecha: str, año: int, mes: int | None) -> bool:
    """
    Verifica si la fecha normalizada cumple el filtro.
    fecha siempre llega como 'YYYY-MM-DD'.
    """
    partes = fecha.split('-')
    año_seq, mes_seq = int(partes[0]), int(partes[1])  # esto cambia "04" → 4
    if año_seq != año:
        return False
    if mes is None:  # Todos los meses
        return True
    return mes_seq == mes

def procesar_archivo(archivo_entrada: str, año: int, mes: int | None) -> tuple[dict, dict]:
    """Filtra secuencias por fecha y separa por presencia de ambigüedades (N)."""
    secuencias_filtradas = {}
    secuencias_con_N = {}
    etiqueta_actual, fecha_actual, secuencia_actual = None, None, []

    with open(archivo_entrada, 'r') as archivo:
        for linea in archivo:
            linea = linea.strip()
            if not linea:
                continue

            if linea.startswith('>'):  # Cabecera FASTA
                # Guardar secuencia previa
                if etiqueta_actual and fecha_actual:
                    secuencia_completa = ''.join(secuencia_actual)
                    if validar_fecha(fecha_actual, año, mes):
                        if 'N' in secuencia_completa:
                            secuencias_con_N[etiqueta_actual] = secuencia_completa
                        else:
                            secuencias_filtradas[etiqueta_actual] = secuencia_completa

                # Nueva cabecera
                etiqueta_actual = linea.split()[0]
                fecha_actual = extraer_fecha(linea)
                secuencia_actual = []

            else:
                if etiqueta_actual:
                    secuencia_actual.append(linea)

        # Procesar última secuencia
        if etiqueta_actual and fecha_actual:
            secuencia_completa = ''.join(secuencia_actual)
            if validar_fecha(fecha_actual, año, mes):
                if 'N' in secuencia_completa:
                    secuencias_con_N[etiqueta_actual] = secuencia_completa
                else:
                    secuencias_filtradas[etiqueta_actual] = secuencia_completa

    return secuencias_filtradas, secuencias_con_N

def guardar_secuencias(secuencias: dict, archivo_salida: str):
    """Guarda secuencias en formato FASTA."""
    with open(archivo_salida, "w") as file:
        for etiqueta, secuencia in secuencias.items():
            file.write(f"{etiqueta}\n{secuencia}\n")

if __name__ == "__main__":
    try:
        mes_str = f", Mes: {MES}" if MES else ", Todos los meses"
        print(f"🔹 Filtrando secuencias de {ARCHIVO_ENTRADA} (Año: {AÑO}{mes_str})...")
        secuencias, secuencias_N = procesar_archivo(ARCHIVO_ENTRADA, AÑO, MES)
        
        guardar_secuencias(secuencias, ARCHIVO_SALIDA)
        print(f"✅ Secuencias sin N guardadas en '{ARCHIVO_SALIDA}'. {len(secuencias)} secuencias")
        
        guardar_secuencias(secuencias_N, ARCHIVO_SALIDA_N)
        print(f"✅ Secuencias con N guardadas en '{ARCHIVO_SALIDA_N}'. {len(secuencias_N)} secuencias")
        
        # Mensaje resumen
        print(f"📊 Total de secuencias procesadas: {len(secuencias) + len(secuencias_N)}")
        
    except FileNotFoundError:
        print(f"❌ Error: No se encontró el archivo '{ARCHIVO_ENTRADA}' o 'parametros.json'.")
    except json.JSONDecodeError:
        print("❌ Error: 'parametros.json' tiene un formato inválido.")
    except KeyError as e:
        print(f"❌ Error: Falta la clave {str(e)} en 'parametros.json'.")
    except Exception as e:
        print(f"❌ Error inesperado: {str(e)}")