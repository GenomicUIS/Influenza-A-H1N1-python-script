from reportlab.lib.pagesizes import letter
from reportlab.lib import colors
from reportlab.pdfgen import canvas
import json
from pathlib import Path

# Cargar configuración desde JSON
with open("parametros.json", "r", encoding="utf-8") as f:
    config = json.load(f)

# Parámetros desde JSON
CEBADORES_FILE = config["cebador"]["conjunto_cebadores"]  # Archivo con los sets de cebadores
COLORS = {
    "directo": getattr(colors, config["pdf"]["color_directo"]),
    "sonda": getattr(colors, config["pdf"]["color_sonda"]),
    "reverso": getattr(colors, config["pdf"]["color_reverso"])
}
OUTPUT_PDF = config["pdf"]["archivo_salida"]
PERIODO = config["filtro"]["periodo"]
MES = config["filtro"]["mes"] if config["filtro"]["mes"] is not None else "todos los meses"

# Selección de consenso según parámetros
usar_consenso = config.get("reporte", {}).get("usar_consenso", "ugene")
if usar_consenso == "biopython":
    CONSENSO_FILE = config["biopython_consensus"]["archivo_salida"]
else:
    CONSENSO_FILE = config["ugene"]["archivo_salida"]

# Diccionario IUPAC para bases ambiguas
iupac_codes = {
    'A': {'A'}, 'T': {'T'}, 'C': {'C'}, 'G': {'G'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'}, 'W': {'A', 'T'},
    'K': {'G', 'T'}, 'M': {'A', 'C'}, 'B': {'C', 'G', 'T'},
    'D': {'A', 'G', 'T'}, 'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'},
    'N': {'A', 'C', 'G', 'T'}
}

iupac_scores = {
    'A': 1.0, 'T': 1.0, 'C': 1.0, 'G': 1.0,
    'R': 0.8, 'Y': 0.8, 'S': 0.8, 'W': 0.8,
    'K': 0.8, 'M': 0.8, 'B': 0.6, 'D': 0.6,
    'H': 0.6, 'V': 0.6, 'N': 0.4
}

def calcular_puntaje_coincidencia(base_consensus, base_primer):
    """Calcula el puntaje de matching entre dos bases considerando IUPAC"""
    if base_consensus == base_primer:
        return 1.0
    
    if base_consensus in iupac_codes and base_primer in iupac_codes:
        interseccion = iupac_codes[base_consensus] & iupac_codes[base_primer]
        if interseccion:
            return min(iupac_scores[base_consensus], iupac_scores[base_primer])
    
    elif base_primer in iupac_codes:
        if base_consensus in iupac_codes[base_primer]:
            return iupac_scores[base_primer]
    
    elif base_consensus in iupac_codes:
        if base_primer in iupac_codes[base_consensus]:
            return iupac_scores[base_consensus]
    
    return 0.0

def read_cebador_sets(file_path):
    """
    Lee sets de cebadores manejando nombres multilínea correctamente.
    Formato esperado:
    19) Nombre muy largo que puede
        continuar en múltiples líneas
    >cebador_directo
    >cebador_sonda
    >cebador_reverso
    """
    sets = []
    current_name = ""
    current_seqs = []
    
    with open(file_path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
            
        # Detectar inicio de nuevo set (línea con número o paréntesis)
        if line[0].isdigit() or line.startswith(')') or (not line.startswith('>') and not current_name):
            if current_name and current_seqs:
                sets.append([current_name] + current_seqs)
            
            current_name = line
            current_seqs = []
            i += 1
            
            # Leer líneas adicionales del nombre (si existen)
            while i < len(lines) and not lines[i].strip().startswith('>') and lines[i].strip():
                current_name += " " + lines[i].strip()
                i += 1
                
        # Leer secuencias de cebadores
        elif line.startswith('>'):
            seq = line[1:].strip().upper()
            if seq:
                current_seqs.append(seq)
            i += 1
            
            # Leer continuación de secuencia multilínea (si existe)
            while i < len(lines) and not lines[i].strip().startswith('>') and lines[i].strip() and not lines[i].strip()[0].isdigit():
                # Continuación de secuencia (sin '>')
                current_seqs[-1] += lines[i].strip().upper()
                i += 1
        else:
            i += 1
    
    # Añadir el último set
    if current_name and current_seqs:
        sets.append([current_name] + current_seqs)
    
    return sets


def find_best_match_for_set(consensus_seq, primer_set):
    """Encuentra mejor coincidencia para un set completo de cebadores"""
    best_match = {
        'set_name': primer_set[0],
        'directo': None,
        'sonda': None,
        'reverso': None,
        'puntaje_total': 0.0,
        'posiciones': None,
        'espaciamiento': None
    }
    
    if len(primer_set) < 3:
        return best_match
    
    directo = primer_set[1]
    sonda = primer_set[2] if len(primer_set) > 2 and primer_set[2] else None
    reverso = primer_set[3] if len(primer_set) > 3 else primer_set[2] if len(primer_set) == 3 and not sonda else None
    
    for i in range(len(consensus_seq) - len(directo) + 1):
        directo_score = sum(calcular_puntaje_coincidencia(consensus_seq[i + j], directo[j]) 
                      for j in range(len(directo)))
        directo_identity = directo_score / len(directo)
        
        if reverso:
            min_reverso_pos = i + len(directo) + 50
            max_reverso_pos = i + len(directo) + 300
            
            for k in range(min_reverso_pos, min(len(consensus_seq) - len(reverso) + 1, max_reverso_pos)):
                reverso_score = sum(calcular_puntaje_coincidencia(consensus_seq[k + l], reverso[l]) 
                              for l in range(len(reverso)))
                reverso_identity = reverso_score / len(reverso)
                
                # Inicializar variables de sonda
                sonda_score = 0.0
                sonda_identity = 0.0  # Inicializar aquí
                sonda_pos = None
                
                if sonda:
                    sonda_pos = i + len(directo) + (k - (i + len(directo))) // 2
                    if sonda_pos + len(sonda) <= k:
                        sonda_score = sum(calcular_puntaje_coincidencia(consensus_seq[sonda_pos + m], sonda[m]) 
                                  for m in range(len(sonda)))
                        sonda_identity = sonda_score / len(sonda)  # Actualizar solo si hay sonda válida
                
                # CALCULAR PUNTAJE TOTAL COMO PROMEDIO
                total_elements = 2  # Directo y reverso siempre existen
                total_score = directo_identity + reverso_identity
                
                if sonda:  # Si hay sonda (aunque su identidad pueda ser 0)
                    total_score += sonda_identity
                    total_elements += 1
                
                average_score = total_score / total_elements
                
                if average_score > best_match['puntaje_total']:
                    best_match = {
                        'set_name': primer_set[0],
                        'directo': {
                            'cebador': directo,
                            'puntaje': directo_identity,
                            'posiciones': i,
                            'puntaje_total': directo_score
                        },
                        'sonda': {
                            'cebador': sonda,
                            'puntaje': sonda_identity,
                            'posiciones': sonda_pos,
                            'puntaje_total': sonda_score
                        } if sonda else None,
                        'reverso': {
                            'cebador': reverso,
                            'puntaje': reverso_identity,
                            'posiciones': k,
                            'puntaje_total': reverso_score
                        },
                        'puntaje_total': average_score,
                        'posiciones': (i, k),
                        'espaciamiento': k - i - len(directo)
                    }
    
    return best_match

def get_comparison_symbol(base_consensus, base_primer):
    """Devuelve un símbolo basado en el nivel de coincidencia manteniendo el puntaje IUPAC original"""
    if base_consensus == base_primer:
        return '|'  # Coincidencia exacta
    
    # Cálculo el puntaje
    score = calcular_puntaje_coincidencia(base_consensus, base_primer)
    
    # Definición de umbrales basados en el puntaje IUPAC
    if score >= 0.8:    # Para bases con puntaje 0.8 o superior (R,Y,S,W,K,M)
        return ':'  
    elif score >= 0.6:  # Para bases con puntaje 0.6 (B,D,H,V)
        return '.'  
    elif score > 0:     # Para bases con puntaje >0 pero <0.6 (solo N con 0.4)
        return '~'  
    else:
        return ' '      # Sin coincidencia

def export_to_pdf(consensus_seq, best_set, archivo_salida=OUTPUT_PDF):
    """Genera un PDF con márgenes de 1 cm y alineación mejorada"""
    c = canvas.Canvas(archivo_salida, pagesize=letter)
    width, height = letter
    
    # Definir márgenes
    MARGIN_LEFT = 30.00
    MARGIN_RIGHT = width - 30.00
    MARGIN_TOP = height - 30.00
    MARGIN_BOTTOM = 30.00
    
    CHAR_WIDTH = 6.0

    # Configuración de estilos
    styles = {
        'title': ('Helvetica-Bold', 14),
        'set_name': ('Helvetica-Bold', 12),
        'primer_name': ('Helvetica-Bold', 10),
        'primer_seq': ('Courier', 10),
        'stats': ('Helvetica', 10),
        'alignment': ('Courier', 10),
        'consensus_title': ('Helvetica-Bold', 12),
        'consensus_seq': ('Courier', 10),
        'pos_num': ('Helvetica', 8),
        'match_title': ('Helvetica-Bold', 10),
        'match_text': ('Courier', 10)
    }
    
    # Encabezado con márgenes
    c.setFont(*styles['title'])
    c.drawString(MARGIN_LEFT, MARGIN_TOP, f"Análisis de Cebadores - Período: {PERIODO}, Mes: {MES}")
    
    # Resultados del mejor set
    y_pos = MARGIN_TOP - 40
    
    if best_set['puntaje_total'] > 0:
        # ==============================================================
        # NUEVO: MANEJO DE NOMBRES LARGOS CON SALTOS DE LÍNEA
        # ==============================================================
        set_name = best_set['set_name']
        max_chars_per_line = 80  # Ajustar según necesidad
        
        # Dividir el nombre en múltiples líneas si es muy largo
        name_lines = []
        while len(set_name) > max_chars_per_line:
            # Encontrar el último espacio dentro del límite
            break_point = set_name.rfind(' ', 0, max_chars_per_line)
            if break_point == -1:
                break_point = max_chars_per_line
            name_lines.append(set_name[:break_point])
            set_name = set_name[break_point:].strip()
        name_lines.append(set_name)
        
        # Dibujar nombre del set en múltiples líneas
        c.setFont(*styles['set_name'])
        for i, line in enumerate(name_lines):
            if i == 0:
                c.drawString(MARGIN_LEFT, y_pos, f"Mejor Set: {line}")
            else:
                c.drawString(MARGIN_LEFT + 80, y_pos, line)  # Sangría para líneas adicionales
            y_pos -= 15
        
        # Puntaje total (en la misma posición que antes)
        c.drawString(MARGIN_RIGHT - 120, MARGIN_TOP, f"Puntaje total: {best_set['puntaje_total']:.2f}")
        y_pos -= 10
        
        # Espaciamiento
        c.setFont(*styles['stats'])
        c.drawString(MARGIN_LEFT, y_pos, f"Espaciamiento directo-reverso: {best_set['espaciamiento']} bases")
        y_pos -= 30
        
        # ==============================================================
        # RESTANTE DEL CÓDIGO ORIGINAL (sin cambios)
        # ==============================================================
        # Mostrar cada cebador del set (directo, sonda, reverso)
        for key in ['directo', 'sonda', 'reverso']:
            data = best_set.get(key)
            if not data:  # 🔹 Saltar si no existe (sonda ausente)
                continue

            # Nombre del cebador
            c.setFillColor(COLORS[key])
            c.setFont(*styles['primer_name'])
            c.drawString(MARGIN_LEFT, y_pos, f"{key.capitalize()}:")
            
            # Secuencia del cebador
            c.setFillColor(colors.black)
            c.setFont(*styles['primer_seq'])
            c.drawString(MARGIN_LEFT + 70, y_pos, data['cebador'])
            y_pos -= 15
            
            # Estadísticas
            c.setFont(*styles['stats'])
            stats = (f"Posición: {data['posiciones']+1}-{data['posiciones']+len(data['cebador'])} | "
                     f"Identidad: {data['puntaje']:.2%} | "
                     f"Puntaje: {data['puntaje_total']:.2f}")
            c.drawString(MARGIN_LEFT, y_pos, stats)
            y_pos -= 20
            
            # COMPARACIÓN DETALLADA
            c.setFont(*styles['match_title'])
            c.drawString(MARGIN_LEFT, y_pos, "Comparación con secuencia consenso:")
            y_pos -= 15
            
            # Obtener segmento del consenso
            start = data['posiciones']
            end = start + len(data['cebador'])
            consensus_segment = consensus_seq[start:end]
            primer_seq = data['cebador']
            
            # Mostrar comparación línea por línea
            line_length = 60
            for i in range(0, len(primer_seq), line_length):
                if y_pos < MARGIN_BOTTOM + 100:
                    c.showPage()
                    y_pos = MARGIN_TOP - 50
                    c.setFont(*styles['match_text'])
                
                current_consensus = consensus_segment[i:i+line_length]
                current_primer = primer_seq[i:i+line_length]
                            
                # Línea de consenso
                c.setFillColor(colors.black)
                c.setFont(*styles['match_text'])
                c.drawString(MARGIN_LEFT + 5, y_pos, "Consenso: ")
                c.drawString(MARGIN_LEFT + 5 + len("Consenso: ") * CHAR_WIDTH, y_pos, current_consensus)
                y_pos -= 12
                
                # Linea de comparacion
                match_symbols = []
                for c_base, p_base in zip(current_consensus, current_primer):
                    match_symbols.append(get_comparison_symbol(c_base, p_base))
                
                for pos, symbol in enumerate(match_symbols):
                    x_pos = MARGIN_LEFT + 5 + len("Consenso: ") * CHAR_WIDTH + pos * CHAR_WIDTH
                    c.drawString(x_pos, y_pos, symbol)
                y_pos -= 12
                
                # Línea de cebador alineada
                c.setFillColor(COLORS[key])
                c.drawString(MARGIN_LEFT, y_pos, f"{key.capitalize():>9}: ")
                c.drawString(MARGIN_LEFT + len(f"{key.capitalize():>9}: ") * CHAR_WIDTH, y_pos, current_primer)
                c.setFillColor(colors.black)
                y_pos -= 15
            
    # Secuencia consenso completa
    c.setFont(*styles['consensus_title'])
    c.drawString(MARGIN_LEFT, y_pos - 20, "Secuencia Consenso Completa:")
    c.setFont(*styles['consensus_seq'])
    
    chars_per_line = 80
    line_height = 12
    margin_y = y_pos - 40
    
    for i in range(0, len(consensus_seq), chars_per_line):
        line = consensus_seq[i:i + chars_per_line]
        y = margin_y - ((i // chars_per_line) * line_height)
        
        if y < MARGIN_BOTTOM:
            c.showPage()
            y = MARGIN_TOP - 50
            c.setFont(*styles['consensus_seq'])
        
        # Números de posición
        c.setFillColor(colors.gray)
        c.setFont(*styles['pos_num'])
        c.drawString(MARGIN_LEFT - 20, y + 3, str(i + 1))
        c.setFillColor(colors.black)
        
        # Bases con resaltamiento
        for j, char in enumerate(line):
            global_pos = i + j
            for key in ['directo', 'sonda', 'reverso']:
                if best_set.get(key):
                    data = best_set[key]
                    if data and data['posiciones'] <= global_pos < data['posiciones'] + len(data['cebador']):
                        c.setFillColor(COLORS[key])
                        c.rect(MARGIN_LEFT + j * 7, y - 2, 7, 10, fill=True, stroke=False)
            
            c.setFillColor(colors.black)
            c.drawString(MARGIN_LEFT + j * 7, y, char)
    
    c.save()
    print(f"\nPDF generado: {archivo_salida}")

def main():
    # === Leer archivo de consenso seleccionado ===
    try:
        with open(CONSENSO_FILE, "r") as f:
            lines = f.readlines()
            consensus_seq = "".join(lines[1:]).replace("\n", "").replace(" ", "") if len(lines) > 1 else ""
    except FileNotFoundError:
        print(f"Error: No se encontró el archivo de consenso {CONSENSO_FILE}")
        return
    
    try:
        sets = read_cebador_sets(CEBADORES_FILE)
    except FileNotFoundError:
        print(f"Error: No se encontró el archivo {CEBADORES_FILE}")
        return
    
    best_set = {
        'set_name': "Ninguno",
        'puntaje_total': 0.0
    }
    
    for primer_set in sets:
        if len(primer_set) >= 3:
            current_set = find_best_match_for_set(consensus_seq, primer_set)
            if current_set['puntaje_total'] > best_set['puntaje_total']:
                best_set = current_set
    
    print("="*70)
    print("RESULTADOS DEL ANÁLISIS DE CEBADORES".center(70))
    print("="*70)
    
    if best_set['puntaje_total'] > 0:
        print(f"\n🔹 MEJOR SET: {best_set['set_name']} (Puntaje total: {best_set['puntaje_total']:.2f})")
        print(f"  Espaciamiento directo-reverso: {best_set['espaciamiento']} bases")
        
        for key in ['directo', 'sonda', 'reverso']:
            if best_set.get(key):
                data = best_set[key]
                print(f"\n  {key.upper()}:")
                print(f"    Cebador: {data['cebador']}")
                print(f"    Posición: {data['posiciones']}")
                print(f"    Identidad: {data['puntaje']:.2%}")
                print(f"    Puntaje: {data['puntaje_total']:.2f}")
    else:
        print("\nNo se encontró ningún set de cebadores con un match adecuado.")
    
    export_to_pdf(consensus_seq, best_set)

if __name__ == "__main__":
    main()