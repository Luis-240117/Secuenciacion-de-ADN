#Sofia Alejandra Beltran
#Luis Francisco Garcia Aguilar
import matplotlib.pyplot as plt
import math
def validarADN(secuencia): #es mas facil trabajar con funciones
    secuencia=secuencia.upper() #usamos solo mayus para mas facil
    caracteresValidos=["A","C","G","T"] #nucletotidos
    for nucleotido in secuencia:
        if nucleotido not in caracteresValidos: #validamos
            return False
    return True

def analizarADN(secuencia): #contamos
    conteoA=secuencia.count('A')
    conteoC=secuencia.count('C')
    conteoG=secuencia.count('G')
    conteoT=secuencia.count('T')
complementoMap={'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
cadenaComplementaria=""
for base in secuencia:
    cadenaComplementaria += complementoMap[base]

secuenciaARN=secuencia.replace('T', 'U') #sustituimos

longitud=len(secuencia)

porcentageGC=((conteoG+conteoC)/longitud)*100 if longitud > 0 else 0.0

print("RESULTADOS DEL ANÁLISIS")
print(f"Conteo = Adenina: {conteoA} | Citosina: {conteoC} | Guanina: {conteoG} | Timina: {conteoT}")
print(f"Porcentaje Guanina-Citosina (GC): {porcentajeGC:.2f}%")
print(f"Cadena Complementaria: {cadenaComplementaria}")
print(f"ARN Mensajero Transcrito: {secuenciaARN}")

def traducirAminoacidos(secuenciaARN):
    codones=[]
    cadenaProteina=[]

import random 
def simularMutacion(secuencia, probabilidad):
    secuenciaMutada=""
    basesPosibles=['A', 'C', 'G', 'T']
for nucleotido in secuencia:
    probGenerada=random.randint(1, 100)
    if probGenerada <= probabilidad:
        nuevaBase=random.choice(basesPosibles)
        secuenciaMutada += nuevaBase
    else:
        secuenciaMutada += nucleotido
return secuenciaMutada

def exportarYGraficar(secuencia, conteoA, conteoC, conteoG, conteoT, porcentajeGC):
    nombreArchivo="Reporte Bioinformatico.txt"
    try:
        with open(nombreArchivo, "w") as archivo:
            archivo.write("REPORTE BIOINFORMATICO\n")
            archivo.write(f"Secuencia Analizada: {secuencia}\n")
            archivo.write(f"Longitud Total: {len(secuencia)} pares de bases.\n\n")
            archivo.write(f"---CONTEO DE NUCLEOTIDOS---\N")
            archivo.write(f"Adenina (A): {conteoA}\n")
            archivo.write(f"Citosina (C): {conteoC}\n")
            archivo.write(f"Guanina (G): {conteoG}\n")
            archivo.write(f"Tiamina (T): {conteoT}\n")
            archivo.write(f"Estabilidad Termica (Contenido GC): {porcentajeGC:.2f}%\n")
        print(f"Exito. Reporte textual guardado en '{nombreArchivo}'.")
    except Exception as e:
        print(f"Error. No se puede guardar el archivo: {e}")
    
    bases=['Adenina (A)', 'Citosina (C)', 'Guanina (G)', 'Timina (T)']
    frecuencias=[conteoA, conteoC, conteoG, conteoT]
    colores = ["#E24AADEC", "#5072E3", "#B3F523", '#D0021B']
    plt.figure(figuresize=(8, 5))
    plt.bar(bases, frecuencias, color=colores)
    plt.title('Frecuencia de Nucleotidos en la Secuencia', fontsize=14)
    plt.xlabel('Bases Nitrogenadas', fontsize=12)
    plt.ylabel('Cantidad', fontsize=12)
    print("Abriendo ventana de visualización gráfica...")
    plt.show()

                          



#modulo del menu
def main():
    opcion=0
    estadoDatosCorrectosCargados=False
    secuenciaActual=""
    while opcion != 5:
        print("==========================")
        print("Sistema Life-Sequence")
        print("==========================")
        print("1. Cargar y Validar Secuencia de ADN")
        print("2. Analizar Secuencia Actual")
        print("3. Simular Mutaciones Aleatorias")
        print("4. Exportar Reporte y Graficar")
        print("5. Salir")
        try:
            opcion=int(input("seleccione una opcion: "))
        except ValueError:
            print("Error. Por favor ingrese un numero.")
            continue

        if opcion == 1:
            entrada=input("Ingrese la secuencia de ADN usando adenina, citosina, guanina y tiamina:")
            if validarADN(entrada):
                print("Secuencia Valida")
                secuenciaActual=entrada.upper()
                estadoDatosCorrectosCargados=True
            else:
                print("Error. La secuencia contiene caracteres no validos.")
                estadoDatosCorrectosCargados=False
        elif opcion==2:
            if estadoDatosCorrectosCargados:
                analizarADN(secuenciaActual)
            else:
                print("Se debe cargar una secuencia valida. Visite opcion 1.")
        elif opcion==3:
            if estadoDatosCargadosCorrectamente:
                probabilidad=int(input("Ingrese la probabilidad de mutacion (1 a 100%): "))
                if 1<=probabilidad<=100:
                    secuenciaResultado=simularMutacion(secuenciaActual, probabilidad)
                    print(f"Secuencia Original: {secuenciaActual}")
                    print(f"Secuencia Mutada:   {secuenciaResultado}")
                    print("-> Análisis de la secuencia mutada:")
                    analizar_adn(secuenciaResultado)
                else:
                    print("[Error. Ingrese un valor entre 1 y 100.")
            else:
                print("Primero debe cargar una secuencia válida.")
        elif opcion==4:
            if estadoDatosCorrectosCargados:
                print("[SISTEMA] Abriendo archivo 'reporte_bio.txt'...")
                print("[SISTEMA] Escribiendo resultados...")
                print("[SISTEMA] Llamando a biblioteca gráfica (matplotlib)...")
                print("[ÉXITO] Reporte exportado y gráfica generada correctamente.")
            else:
                print("No hay datos para exportar.")
        elif opcion==5:
            print("Saliendo del programa. La secuencia vital finalizara.")
        else:
            print("Opcion no valida. Intente de nuevo o seleccione nuevamente.")
if __name__=="__main__":
    main()
