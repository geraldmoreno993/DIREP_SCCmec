from Bio import SeqIO

def buscar_repeticiones_directas_con_posiciones(secuencia, min_len=15, max_len=30, min_dist=10000, max_dist=60000):
    repes_filtradas = []

    for length in range(min_len, max_len + 1):
        substrings = {}
        for i in range(len(secuencia) - length + 1):
            subseq = secuencia[i:i + length]
            if subseq not in substrings:
                substrings[subseq] = [i]
            else:
                substrings[subseq].append(i)

        for subseq, posiciones in substrings.items():
            if len(posiciones) < 2 or len(posiciones) > 10:  # descarta repes muy frecuentes
                continue

            for i in range(len(posiciones)):
                for j in range(i + 1, len(posiciones)):
                    dist = abs(posiciones[j] - posiciones[i])
                    if min_dist <= dist <= max_dist:
                        repes_filtradas.append({
                            'secuencia': subseq,
                            'longitud': length,
                            'inicio_1': posiciones[i],
                            'fin_1': posiciones[i] + length,
                            'inicio_2': posiciones[j],
                            'fin_2': posiciones[j] + length,
                            'distancia': dist
                        })

    return repes_filtradas


# === PARTE PRINCIPAL ===
from pathlib import Path
from sys import argv
from Bio import SeqIO

# Usa archivo FASTA pasado por argumento o por defecto
if len(argv) > 1:
    fasta = argv[1]
else:
    fasta = "R10_chop_contig_3_sccmec.fasta"  # cambiar si deseas otro nombre por defecto

# Verifica existencia
if not Path(fasta).exists():
    print(f"[ERROR] No se encontró el archivo: {fasta}")
    exit(1)

# Cargar secuencia
registro = SeqIO.read(fasta, "fasta")
seq = str(registro.seq)

# Buscar repeticiones
resultado = buscar_repeticiones_directas_con_posiciones(seq)

# Mostrar resultados
print(f"{'Secuencia':<30} {'Long':<5} {'Inicio1':<8} {'Fin1':<8} {'Inicio2':<8} {'Fin2':<8} {'Distancia'}")
for r in resultado:
    print(f"{r['secuencia']:<30} {r['longitud']:<5} {r['inicio_1']:<8} {r['fin_1']:<8} {r['inicio_2']:<8} {r['fin_2']:<8} {r['distancia']}")


