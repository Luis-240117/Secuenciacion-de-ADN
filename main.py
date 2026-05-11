#Luis Francisco Garcia Aguilar
#Sofia Alejandra Beltran Reyes

import streamlit as st        #Importa streamlit para poder construir la interfaz
import matplotlib.pyplot as plt  #Importa matplotlib para generar los graficos.
import random                 

# Diccionario local con el cual revisamos cada codon, y a su vez marcamos tambien los inicios y los fines de cada secuencia.
codigoGeneticoCompleto = {
    'UUU': 'Fenilalanina', 'UUC': 'Fenilalanina', 'UUA': 'Leucina', 'UUG': 'Leucina',
    'CUU': 'Leucina', 'CUC': 'Leucina', 'CUA': 'Leucina', 'CUG': 'Leucina',
    'AUU': 'Isoleucina', 'AUC': 'Isoleucina', 'AUA': 'Isoleucina', 'AUG': 'Metionina(Inicio)',  #AUG codon de inicio universal.
    'GUU': 'Valina', 'GUC': 'Valina', 'GUA': 'Valina', 'GUG': 'Valina',
    'UCU': 'Serina', 'UCC': 'Serina', 'UCA': 'Serina', 'UCG': 'Serina', 'AGU': 'Serina', 'AGC': 'Serina',
    'CCU': 'Prolina', 'CCC': 'Prolina', 'CCA': 'Prolina', 'CCG': 'Prolina',
    'ACU': 'Treonina', 'ACC': 'Treonina', 'ACA': 'Treonina', 'ACG': 'Treonina',
    'GCU': 'Alanina', 'GCC': 'Alanina', 'GCA': 'Alanina', 'GCG': 'Alanina',
    'UAU': 'Tirosina', 'UAC': 'Tirosina', 'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop',  #UAA, UAG y UGA codones de parada.
    'CAU': 'Histidina', 'CAC': 'Histidina', 'CAA': 'Glutamina', 'CAG': 'Glutamina',
    'AAU': 'Asparagina', 'AAC': 'Asparagina', 'AAA': 'Lisina', 'AAG': 'Lisina',
    'GAU': 'Ácido Aspártico', 'GAC': 'Ácido Aspártico', 'GAA': 'Ácido Glutámico', 'GAG': 'Ácido Glutámico',
    'UGU': 'Cisteína', 'UGC': 'Cisteína', 'UGG': 'Triptófano',
    'CGU': 'Arginina', 'CGC': 'Arginina', 'CGA': 'Arginina', 'CGG': 'Arginina', 'AGA': 'Arginina', 'AGG': 'Arginina',
    'GGU': 'Glicina', 'GGC': 'Glicina', 'GGA': 'Glicina', 'GGG': 'Glicina'
}

#Verifica que la secuencia ingresada solo contenga nucleotidos validos de ADN
def validarADN(secuencia):
    secuencia = secuencia.upper()           #Convierte toda la cadena a maysculas para tener un estandar en el manejo de la entrada.
    caracteresValidos = ["A", "C", "G", "T"]  #Lista de los cuatro nucleotidos validos.
    for nucleotido in secuencia:            #Recorre cada caracter.
        if nucleotido not in caracteresValidos:
            return False                    
    return True                             

#Realiza el analisis basico de la secuencia.
def analizarADN(secuencia):
    secuencia = secuencia.upper()
    conteoA = secuencia.count('A') #Cuenta cuantas Adeninas tenemos.
    conteoC = secuencia.count('C') #Cuenta cuantas Citosinas tenemos.
    conteoG = secuencia.count('G') #Cuenta cuantas Guaninas tenemos.
    conteoT = secuencia.count('T') #Cuenta cuantas Timinas tenemos.
    
    #Diccionario de emparejamiento de bases segun las reglas de Watson y Crick.
    complementoMap = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    #Genera la cadena complementaria sustituyendo cada base por su par complementario
    cadenaComplementaria = "".join([complementoMap[base] for base in secuencia])
    
    #Transcripcion
    secuenciaARN = secuencia.replace('T', 'U')
    
    longitud = len(secuencia) # Calcula la longitud de pares de bases.
    
    #Calcula el porcentaje de Guanina + Citosina respecto al total.
    porcentajeGC = ((conteoG + conteoC) / longitud) * 100 if longitud > 0 else 0.0 #Evitamos una division sobre cero.
    
    #Retorna todos los valores calculados.
    return conteoA, conteoC, conteoG, conteoT, porcentajeGC, cadenaComplementaria, secuenciaARN

#Calcula propiedades biofisicas adicionales de la secuencia de ADN.
def calcularMetricasExtra(secuencia, conteoA, conteoC, conteoG, conteoT):
    #Temperatura de fusion usando la regla de Wallace.
    #Asigna 2°C por cada par A-T y 4°C por cada par G-C, basado en la cantidad de enlaces que forman entre si cada uno de los pares.
    tm = (2 * (conteoA + conteoT)) + (4 * (conteoG + conteoC))
    
    #Peso molecular aproximado en gramos/mol para ADN de cadena sencilla
    #Se resta 61.96 para corregir la perdida de agua durante la formacion de enlaces fosfodiester.
    pesoMolecular = (conteoA * 313.2) + (conteoC * 289.2) + (conteoG * 329.2) + (conteoT * 304.2) - 61.96
    
    #Reutiliza el mapa de complementos para calcular el complemento inverso
    complementoMap = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    #Genera la cadena complementaria base a base
    cadenaComplementaria = "".join([complementoMap[base] for base in secuencia])
    
    #Invierte la cadena complementaria para obtener la dirección 3' → 5'
    complementoInverso = cadenaComplementaria[::-1]
    
    #Nos retorna los calculos.
    return tm, pesoMolecular, complementoInverso

#Simula el proceso de traduccion, lee el ARN mensajero en tripletes.
def traducirAminoacidos(secuenciaARN):
    codones = [] #Lista que almacenara cada codon leído.
    cadenaProteina = [] #Lista que almacenara el aminoacido correspondiente a cada codon.
    
    #Busca el primer codon de inicio AUG en la secuencia.
    inicioCodon = secuenciaARN.find('AUG')
    
    #Si no se encuentra AUG comenzamos a leer desde el inicio de la secuencia ingresada.
    if inicioCodon == -1:
        inicioCodon = 0
        
    #Recorre la secuencia desde el codón de inicio en pasos de 3.
    #Se detiene 2 posiciones antes del final para garantizar que siempre haya un codon completo.
    for i in range(inicioCodon, len(secuenciaARN) - 2, 3):
        codonActual = secuenciaARN[i:i+3] #Extrae el trio de nucleotidos en la posicion actual.
        codones.append(codonActual) # Agrega el codon a la lista de codones procesados.
        
        #Busca el aminosacido correspondiente al codon en el diccionario genetico.
        aminoacidoResultante = codigoGeneticoCompleto.get(codonActual, 'Desconocido') #Si el codon no esta en el diccionario no regresara un desconocido.
        cadenaProteina.append(aminoacidoResultante) #Agrega el aminoacido a la cadena polipeptidica generada.
        
        #Al encontrarnos con un codon de parada la traduccion se detiene como en las celulas.
        if aminoacidoResultante == 'Stop':
            break
                
    #Nos devuelve los codones leidos y los aminoacidos producidos.
    return codones, cadenaProteina

#Simula mutaciones de tipo sustitucion sobre la secuencia de ADN, y cada nucleotido tiene una probabilidad de ser reemplazado por una base aleatoria.
def simularMutacion(secuencia, probabilidad):
    secuenciaMutada = "" #Cadena vacia que acumulara la secuencia.
    basesPosibles = ['A', 'C', 'G', 'T'] #Cuatro bases probables.
    
    for nucleotido in secuencia: #Recorremos la secuencia nucleotido por nucleotido.
        probGenerada = random.randint(1, 100) #Genera un numero aleatorio entre 1 y 100 para simular la probabilidad.
        
        if probGenerada <= probabilidad: #Si el numero generado cae dentro del umbral de probabilidad sustituye el nucleotido por una base aleatoria
            secuenciaMutada += random.choice(basesPosibles)
        else:
            secuenciaMutada += nucleotido
            
    return secuenciaMutada #Nos regresa la secuencia con las mutaciones aplicadas.

#Formatea una secuencia de ADN al estandar FASTA usado en bioinformatica.
def generarFASTA(secuencia, nombreSecuencia="Secuencia_Analizada_Life_Sequence"):
    formatoFasta = f">{nombreSecuencia}\n"  #Nos aseguramos de comenzar con ">".
    
    #Divide la secuencia en fragmentos de 80 caracteres para cumplir con el formato FASTA.
    fragmentosFasta = [secuencia[i:i+80] for i in range(0, len(secuencia), 80)]
    
    #Une todos los fragmentos con saltos de linea para generar el cuerpo del archivo.
    formatoFasta += "\n".join(fragmentosFasta)
    
    return formatoFasta #Nos regresa el formato.


#Configura el titulo en la pestaña del navegador, ícono y diseño amplio.
st.set_page_config(page_title="Sistema Life-Sequence", page_icon="🧬", layout="wide")

st.title("Sistema Life-Sequence") #Nos permite mostrar el titulo de la pagina en la parte superior.
st.markdown("Plataforma interactiva para el análisis bioinformático de secuencias de ADN.")  #Subtitulo

#Inicializa la variable de sesion que almacena la secuencia actual del usuario.
if "secuenciaActual" not in st.session_state:
    st.session_state.secuenciaActual = "" #Cadena vacia.

#Toda la logica de carga de datos se coloca en la barra lateral para poder tener mas espacio para la informacion.
st.sidebar.header("1. Carga de Datos")  #Encabezado del panel lateral.

#Area de texto para ingresar la secuencia.
entradaAdn = st.sidebar.text_area("Ingrese la secuencia de ADN:", value=st.session_state.secuenciaActual)

#Boton que desencadena la validacion y carga de la secuencia ingresada
if st.sidebar.button("Validar y Cargar Secuencia"):
    if entradaAdn: #Verifica que el campo de texto no esté vacio.
        #Limpia la entrada.
        entradaLimpia = entradaAdn.replace(" ", "").replace("\n", "").upper()
        
        if validarADN(entradaLimpia): #Llama a la función de validacion.
            st.session_state.secuenciaActual = entradaLimpia #Guarda la secuencia valida en la sesion.
            st.sidebar.success("¡Secuencia valida y cargada correctamente!") #Mensaje de exito.
        else:
            #Si la secuencia contiene caracteres invalidos, muestra un mensaje de error en rojo.
            st.sidebar.error("Error: La secuencia contiene caracteres no validos. Solo use A, C, G, T.")
    else:
        #Si el campo esta vacio al presionar el boton, muestra una advertencia en amarillo.
        st.sidebar.warning("Por favor, ingrese una secuencia primero.")

#Solo muestra el contenido principal si ya hay una secuencia cargada y validada en la sesion.
if st.session_state.secuenciaActual:
    #Crea tres pestañas navegables para organizar las funciones del sistema.
    tabUno, tabDos, tabTres = st.tabs(["Analisis Principal", "Traduccion a Proteinas", "Simulacion de Mutaciones"])
    
    secuenciaEvaluada = st.session_state.secuenciaActual   # Referencia local a la secuencia guardada en sesión
    
    #Ejecuta el analisis base una sola vez y comparte los resultados entre todas las pestañas.
    conteoAdenina, conteoCitosina, conteoGuanina, conteoTimina, porcentajeGlobalGC, cadenaComp, arnMensajero = analizarADN(secuenciaEvaluada)
    
    #Calcula las metricas biofisicas.
    temperaturaFusion, pesoMol, complementoInv = calcularMetricasExtra(secuenciaEvaluada, conteoAdenina, conteoCitosina, conteoGuanina, conteoTimina)
    
    #Pestaña Uno
    with tabUno:
        st.header("Resultados del Analisis Estructural")  #Encabezado de la pestaña
        st.write(f"**Secuencia Original:** `{secuenciaEvaluada}`")  #Muestra la secuencia en formato monoespaciado.
        
        #Crea 5 columnas iguales para mostrar los conteos de cada nucleotido y el porcentaje GC.
        columnaUno, columnaDos, columnaTres, columnaCuatro, columnaCinco = st.columns(5)
        columnaUno.metric("Adenina (A)", conteoAdenina) #Total de Adeninas.
        columnaDos.metric("Timina (T)", conteoTimina) #Total de Timinas.
        columnaTres.metric("Citosina (C)", conteoCitosina) #Total de Citosinas.
        columnaCuatro.metric("Guanina (G)", conteoGuanina) #Total de Guaninas.
        columnaCinco.metric("Porcentaje GC", f"{porcentajeGlobalGC:.2f}%")  #Porcentaje GC formateado a 2 decimales
        
        st.markdown("---") #Linea horizontal divisora.
        
        #Seccion de propiedades biofisicas calculadas.
        st.subheader("Propiedades Biofisicas") #Subtitulo para esta nueva seccion.
        
        #Tres columnas para mostrar temperatura de fusion, peso molecular y longitud.
        colMetricaUno, colMetricaDos, colMetricaTres = st.columns(3)
        colMetricaUno.metric("Temperatura de Fusion (Tm)", f"{temperaturaFusion} °C") #Temperatura en grados Celsius.
        colMetricaDos.metric("Peso Molecular Estimado", f"{pesoMol:.2f} g/mol") #Peso con 2 decimales.
        colMetricaTres.metric("Longitud Total", f"{len(secuenciaEvaluada)} pb") #Longitud en pares de bases
        
        st.markdown("---")  #Otra division.
        
        #Muestra las cadenas derivadas de la secuencia original en formato monoespaciado.
        st.write(f"**Cadena Complementaria:** `{cadenaComp}`") #Hebra complementaria.
        st.write(f"**Complemento Inverso (3' a 5'):** `{complementoInv}`") # Complemento leido en sentido antiparalelo.
        st.write(f"**ARN Mensajero Transcrito:** `{arnMensajero}`") #ARN resultante de la transcripcion.
        
        st.markdown("---") #Linea divisoria.
        
        # Divide el area en dos columnas.
        colGrafica, colDescargas = st.columns([2, 1])
        
        with colGrafica:
            st.subheader("Frecuencia de Nucleotidos") #Titulo de la grafica.
            
            #Etiquetas del eje X de la grafica de barras.
            basesNitrogenadas = ['Adenina (A)', 'Citosina (C)', 'Guanina (G)', 'Timina (T)']
            
            #Valores de altura de cada barra segun los conteos obtenidos
            frecuenciasBases = [conteoAdenina, conteoCitosina, conteoGuanina, conteoTimina]
            
            #Colores para cada barra.
            coloresGrafica = ["#E24ADE", "#5072E3", "#B3F523", "#DF9606"]
            
            #Crea la figura y el eje de matplotlib.
            figuraGrafico, ejeGrafico = plt.subplots(figsize=(6, 4))
            
            #Dibuja la grafica de barras con las bases, frecuencias y colores.
            ejeGrafico.bar(basesNitrogenadas, frecuenciasBases, color=coloresGrafica)
            
            ejeGrafico.set_ylabel('Cantidad') #Etiqueta del eje Y.
            st.pyplot(figuraGrafico) #Muestra la figura de matplotlib en la interfaz.
            
        with colDescargas:
            st.subheader("Exportar Datos") #Titulo de la seccion de datos que podemos exportar.
            st.write("Guarda tus resultados para analisis posteriores.") #Instruccion para que los usuarios sepan que hacer.
            
            #Construye el contenido del reporte en texto plano con todos los resultados calculados
            reporteTxt = f"""REPORTE BIOINFORMATICO - LIFE SEQUENCE
======================================
Secuencia Original: {secuenciaEvaluada}
Longitud: {len(secuenciaEvaluada)} pb
Porcentaje GC: {porcentajeGlobalGC:.2f}%
Temperatura de Fusion (Tm): {temperaturaFusion} °C
Peso Molecular: {pesoMol:.2f} g/mol
ARN Mensajero: {arnMensajero}
Complemento Inverso: {complementoInv}
======================================"""
            
            #Boton que permite al usuario descargar el reporte como archivo .txt.
            st.download_button("Descargar Reporte (.txt)", data=reporteTxt, file_name="Reporte_Bioinformatico.txt", mime="text/plain", use_container_width=True)
            
            #Genera el contenido en formato FASTA.
            fastaTxt = generarFASTA(secuenciaEvaluada)
            
            #Boton que permite descargar la secuencia en formato FASTA
            st.download_button("Descargar formato FASTA", data=fastaTxt, file_name="secuencia.fasta", mime="text/plain", use_container_width=True)

    #Pestaña Dos
    with tabDos:
        st.header("Prediccion de Proteinas (ORF)")
        st.write(f"**ARN Mensajero (Plantilla):** `{arnMensajero}`") #Muestra el ARN que se va a traducir.
        
        #Llama a la funcion de traducción y obtiene los valores que retorna la funcion.
        listaCodones, listaProteina = traducirAminoacidos(arnMensajero)
        
        if listaProteina: #Verifica que se haya producido al menos un aminoacido.
            st.success("Traduccion completada con exito.")  #Mensaje de exito en verde.
            
            st.write("**Secuencia de Codones procesados:**")
            #Muestra los codones separados por guiones.
            st.code(" - ".join(listaCodones))
            
            st.write("**Secuencia de Aminoacidos (Cadena Polipeptidica):**")
            #Muestra los aminoacidos.
            st.info(" - ".join(listaProteina))
            
            #Si la cadena termino con un codon de parada muestra el mensaje para que el usuario lo sepa.
            if 'Stop' in listaProteina:
                st.caption("Nota: La traduccion se detuvo al encontrar un codon de parada (Stop).")
        else:
            #Si no se generaron aminoacidos muestra una advertencia.
            st.warning("No se generaron aminoacidos. Asegurese de que la secuencia tenga la longitud suficiente.")

    #Pestaña Tres
    with tabTres:
        st.header("Simulacion de Mutaciones Aleatorias") #Titulo de la pestaña.
        st.markdown("Aplica una tasa de mutacion puntual (sustitucion) a la secuencia actual.") #Descripcion del proceso.
        
        #Control deslizante para que el usuario ajuste la probabilidad de mutacion por nucleotido.
        probabilidadMutacion = st.slider("Probabilidad de mutacion por nucleotido (%)", min_value=1, max_value=100, value=10)
        
        #Boton que ejecuta la simulacion usando la probabilidad seleccionada por el usuario.
        if st.button("Generar Mutacion"):
            #Llama a la funcion y obtiene los valores de retorno.
            secuenciaResultado = simularMutacion(secuenciaEvaluada, probabilidadMutacion)
            
            #Muestra ambas secuencias para compararlas.
            st.write(f"**Secuencia Original:** `{secuenciaEvaluada}`")
            st.write(f"**Secuencia Mutada:**   `{secuenciaResultado}`")
            
            #Cuenta cuantas posiciones son distintas entre la secuencia original y la mutada.
            cantidadDiferencias = sum(1 for a, b in zip(secuenciaEvaluada, secuenciaResultado) if a != b)
            
            #Calcula el porcentaje real de mutaciones ocurridas respecto a la longitud total.
            porcentajeDiferencia = (cantidadDiferencias / len(secuenciaEvaluada)) * 100
            
            #Muestra el resumen de mutaciones en un recuadro de advertencia amarillo.
            st.warning(f"Se detectaron **{cantidadDiferencias}** mutaciones ({porcentajeDiferencia:.1f}%) respecto a la secuencia original.")
else:
    #Si aun no se ha cargado ninguna secuencia, mostramos un mensaje guiar al usuario.
    st.info("👈 Por favor, ingresa y valida una secuencia de ADN valida en el panel lateral para comenzar el anlisis.")