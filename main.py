#Luis Francisco Garcia Aguilar
#Sofia Alejandra Beltran Reyes

import streamlit as st
import matplotlib.pyplot as plt
import random

codigoGeneticoCompleto = {
    'UUU': 'Fenilalanina', 'UUC': 'Fenilalanina', 'UUA': 'Leucina', 'UUG': 'Leucina',
    'CUU': 'Leucina', 'CUC': 'Leucina', 'CUA': 'Leucina', 'CUG': 'Leucina',
    'AUU': 'Isoleucina', 'AUC': 'Isoleucina', 'AUA': 'Isoleucina', 'AUG': 'Metionina(Inicio)',
    'GUU': 'Valina', 'GUC': 'Valina', 'GUA': 'Valina', 'GUG': 'Valina',
    'UCU': 'Serina', 'UCC': 'Serina', 'UCA': 'Serina', 'UCG': 'Serina', 'AGU': 'Serina', 'AGC': 'Serina',
    'CCU': 'Prolina', 'CCC': 'Prolina', 'CCA': 'Prolina', 'CCG': 'Prolina',
    'ACU': 'Treonina', 'ACC': 'Treonina', 'ACA': 'Treonina', 'ACG': 'Treonina',
    'GCU': 'Alanina', 'GCC': 'Alanina', 'GCA': 'Alanina', 'GCG': 'Alanina',
    'UAU': 'Tirosina', 'UAC': 'Tirosina', 'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop',
    'CAU': 'Histidina', 'CAC': 'Histidina', 'CAA': 'Glutamina', 'CAG': 'Glutamina',
    'AAU': 'Asparagina', 'AAC': 'Asparagina', 'AAA': 'Lisina', 'AAG': 'Lisina',
    'GAU': 'Ácido Aspártico', 'GAC': 'Ácido Aspártico', 'GAA': 'Ácido Glutámico', 'GAG': 'Ácido Glutámico',
    'UGU': 'Cisteína', 'UGC': 'Cisteína', 'UGG': 'Triptófano',
    'CGU': 'Arginina', 'CGC': 'Arginina', 'CGA': 'Arginina', 'CGG': 'Arginina', 'AGA': 'Arginina', 'AGG': 'Arginina',
    'GGU': 'Glicina', 'GGC': 'Glicina', 'GGA': 'Glicina', 'GGG': 'Glicina'
}


def validarADN(secuencia):
    secuencia = secuencia.upper()
    caracteresValidos = ["A", "C", "G", "T"]
    for nucleotido in secuencia:
        if nucleotido not in caracteresValidos:
            return False
    return True

def analizarADN(secuencia):
    secuencia = secuencia.upper()
    conteoA = secuencia.count('A')
    conteoC = secuencia.count('C')
    conteoG = secuencia.count('G')
    conteoT = secuencia.count('T')
    
    complementoMap = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    cadenaComplementaria = "".join([complementoMap[base] for base in secuencia])
    secuenciaARN = secuencia.replace('T', 'U')
    longitud = len(secuencia)
    porcentajeGC = ((conteoG + conteoC) / longitud) * 100 if longitud > 0 else 0.0
    
    return conteoA, conteoC, conteoG, conteoT, porcentajeGC, cadenaComplementaria, secuenciaARN

def calcularMetricasExtra(secuencia, conteoA, conteoC, conteoG, conteoT):
    # Temperatura de Fusión (Regla de Wallace)
    tm = (2 * (conteoA + conteoT)) + (4 * (conteoG + conteoC))
    
    # Peso molecular aproximado en g/mol para ADN monocatenario
    pesoMolecular = (conteoA * 313.2) + (conteoC * 289.2) + (conteoG * 329.2) + (conteoT * 304.2) - 61.96
    
    # Complemento Inverso
    complementoMap = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    cadenaComplementaria = "".join([complementoMap[base] for base in secuencia])
    complementoInverso = cadenaComplementaria[::-1] # Invierte el orden de 3' a 5'
    
    return tm, pesoMolecular, complementoInverso

def traducirAminoacidos(secuenciaARN):
    codones = []
    cadenaProteina = []
    
    # Busca un codon de inicio.
    inicioCodon = secuenciaARN.find('AUG')
    if inicioCodon == -1:
        inicioCodon = 0 
        
    for i in range(inicioCodon, len(secuenciaARN) - 2, 3):
        codonActual = secuenciaARN[i:i+3]
        codones.append(codonActual)
        aminoacidoResultante = codigoGeneticoCompleto.get(codonActual, 'Desconocido')
        cadenaProteina.append(aminoacidoResultante)
        
        if aminoacidoResultante == 'Stop':
            break # Termina al encontrar un codon de parada.
                
    return codones, cadenaProteina

def simularMutacion(secuencia, probabilidad):
    secuenciaMutada = ""
    basesPosibles = ['A', 'C', 'G', 'T']
    for nucleotido in secuencia:
        probGenerada = random.randint(1, 100)
        if probGenerada <= probabilidad:
            secuenciaMutada += random.choice(basesPosibles)
        else:
            secuenciaMutada += nucleotido
    return secuenciaMutada

def generarFASTA(secuencia, nombreSecuencia="Secuencia_Analizada_Life_Sequence"):
    formatoFasta = f">{nombreSecuencia}\n"
    fragmentosFasta = [secuencia[i:i+80] for i in range(0, len(secuencia), 80)]
    formatoFasta += "\n".join(fragmentosFasta)
    return formatoFasta

# ==========================================
# 3. INTERFAZ GRÁFICA CON STREAMLIT
# ==========================================
st.set_page_config(page_title="Sistema Life-Sequence", page_icon="🧬", layout="wide")

st.title("🧬 Sistema Life-Sequence")
st.markdown("Plataforma interactiva para el análisis bioinformático de secuencias de ADN.")

# Inicializar variable de sesion
if "secuenciaActual" not in st.session_state:
    st.session_state.secuenciaActual = ""

# --- PANEL LATERAL ---
st.sidebar.header("1. Carga de Datos")
entradaAdn = st.sidebar.text_area("Ingrese la secuencia de ADN:", value=st.session_state.secuenciaActual)

if st.sidebar.button("Validar y Cargar Secuencia"):
    if entradaAdn:
        entradaLimpia = entradaAdn.replace(" ", "").replace("\n", "").upper()
        if validarADN(entradaLimpia):
            st.session_state.secuenciaActual = entradaLimpia
            st.sidebar.success("¡Secuencia válida y cargada correctamente!")
        else:
            st.sidebar.error("Error: La secuencia contiene caracteres no válidos. Solo use A, C, G, T.")
    else:
        st.sidebar.warning("Por favor, ingrese una secuencia primero.")

# --- PESTAÑAS PRINCIPALES ---
if st.session_state.secuenciaActual:
    tabUno, tabDos, tabTres = st.tabs(["📊 Análisis Principal", "🧪 Traducción a Proteínas", "🧬 Simulación de Mutaciones"])
    
    secuenciaEvaluada = st.session_state.secuenciaActual
    
    # Análisis de los datos base
    conteoAdenina, conteoCitosina, conteoGuanina, conteoTimina, porcentajeGlobalGC, cadenaComp, arnMensajero = analizarADN(secuenciaEvaluada)
    temperaturaFusion, pesoMol, complementoInv = calcularMetricasExtra(secuenciaEvaluada, conteoAdenina, conteoCitosina, conteoGuanina, conteoTimina)
    
    # ----------------------------------------------------------------------
    # PESTAÑA 1: ANÁLISIS PRINCIPAL
    # ----------------------------------------------------------------------
    with tabUno:
        st.header("Resultados del Análisis Estructural")
        st.write(f"**Secuencia Original:** `{secuenciaEvaluada}`")
        
        #Metricas nucleotidos.
        columnaUno, columnaDos, columnaTres, columnaCuatro, columnaCinco = st.columns(5)
        columnaUno.metric("Adenina (A)", conteoAdenina)
        columnaDos.metric("Timina (T)", conteoTimina)
        columnaTres.metric("Citosina (C)", conteoCitosina)
        columnaCuatro.metric("Guanina (G)", conteoGuanina)
        columnaCinco.metric("Porcentaje GC", f"{porcentajeGlobalGC:.2f}%")
        
        st.markdown("---")
        
        #Metricas Fisicoquimicas
        st.subheader("Propiedades Biofísicas")
        colMetricaUno, colMetricaDos, colMetricaTres = st.columns(3)
        colMetricaUno.metric("Temperatura de Fusión (Tm)", f"{temperaturaFusion} °C")
        colMetricaDos.metric("Peso Molecular Estimado", f"{pesoMol:.2f} g/mol")
        colMetricaTres.metric("Longitud Total", f"{len(secuenciaEvaluada)} pb")
        
        st.markdown("---")
        
        # Textos de cadenas relacionadas.
        st.write(f"**Cadena Complementaria:** `{cadenaComp}`")
        st.write(f"**Complemento Inverso (3' a 5'):** `{complementoInv}`")
        st.write(f"**ARN Mensajero Transcrito:** `{arnMensajero}`")
        
        st.markdown("---")
        
        # Graficas y botones de descarga.
        colGrafica, colDescargas = st.columns([2, 1])
        
        with colGrafica:
            st.subheader("Frecuencia de Nucleótidos")
            basesNitrogenadas = ['Adenina (A)', 'Citosina (C)', 'Guanina (G)', 'Timina (T)']
            frecuenciasBases = [conteoAdenina, conteoCitosina, conteoGuanina, conteoTimina]
            coloresGrafica = ["#E24ADE", "#5072E3", "#B3F523", '#D0021B']
            
            figuraGrafico, ejeGrafico = plt.subplots(figsize=(6, 4))
            ejeGrafico.bar(basesNitrogenadas, frecuenciasBases, color=coloresGrafica)
            ejeGrafico.set_ylabel('Cantidad')
            st.pyplot(figuraGrafico)
            
        with colDescargas:
            st.subheader("Exportar Datos")
            st.write("Guarda tus resultados para análisis posteriores.")
            
            # Generar TXT
            reporteTxt = f"""REPORTE BIOINFORMATICO - LIFE SEQUENCE
======================================
Secuencia Original: {secuenciaEvaluada}
Longitud: {len(secuenciaEvaluada)} pb
Porcentaje GC: {porcentajeGlobalGC:.2f}%
Temperatura de Fusión (Tm): {temperaturaFusion} °C
Peso Molecular: {pesoMol:.2f} g/mol
ARN Mensajero: {arnMensajero}
Complemento Inverso: {complementoInv}
======================================"""
            
            st.download_button("📥 Descargar Reporte (.txt)", data=reporteTxt, file_name="Reporte_Bioinformatico.txt", mime="text/plain", use_container_width=True)
            
            # Generar FASTA
            fastaTxt = generarFASTA(secuenciaEvaluada)
            st.download_button("🧬 Descargar formato FASTA", data=fastaTxt, file_name="secuencia.fasta", mime="text/plain", use_container_width=True)

    # ----------------------------------------------------------------------
    # PESTAÑA 2: TRADUCCIÓN A PROTEÍNAS
    # ----------------------------------------------------------------------
    with tabDos:
        st.header("Predicción de Proteínas (ORF)")
        st.write(f"**ARN Mensajero (Plantilla):** `{arnMensajero}`")
        
        listaCodones, listaProteina = traducirAminoacidos(arnMensajero)
        
        if listaProteina:
            st.success("Traducción completada con éxito.")
            st.write("**Secuencia de Codones procesados:**")
            st.code(" - ".join(listaCodones))
            
            st.write("**Secuencia de Aminoácidos (Cadena Polipeptídica):**")
            st.info(" - ".join(listaProteina))
            
            if 'Stop' in listaProteina:
                st.caption("Nota: La traducción se detuvo al encontrar un codón de parada (Stop).")
        else:
            st.warning("No se generaron aminoácidos. Asegúrese de que la secuencia tenga la longitud suficiente.")

    # ----------------------------------------------------------------------
    # PESTAÑA 3: MUTACIONES
    # ----------------------------------------------------------------------
    with tabTres:
        st.header("Simulación de Mutaciones Aleatorias")
        st.markdown("Aplica una tasa de mutación puntual (sustitución) a la secuencia actual.")
        
        probabilidadMutacion = st.slider("Probabilidad de mutación por nucleótido (%)", min_value=1, max_value=100, value=10)
        
        if st.button("Generar Mutación"):
            secuenciaResultado = simularMutacion(secuenciaEvaluada, probabilidadMutacion)
            
            st.write(f"**Secuencia Original:** `{secuenciaEvaluada}`")
            st.write(f"**Secuencia Mutada:**   `{secuenciaResultado}`")
            
            # Contar diferencias reales generadas
            cantidadDiferencias = sum(1 for a, b in zip(secuenciaEvaluada, secuenciaResultado) if a != b)
            porcentajeDiferencia = (cantidadDiferencias / len(secuenciaEvaluada)) * 100
            
            st.warning(f"Se detectaron **{cantidadDiferencias}** mutaciones ({porcentajeDiferencia:.1f}%) respecto a la secuencia original.")
else:
    st.info("👈 Por favor, ingresa y valida una secuencia de ADN válida en el panel lateral para comenzar el análisis.")