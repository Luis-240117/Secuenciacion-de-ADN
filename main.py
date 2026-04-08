# Sofia Alejandra Beltran
# Luis Francisco Garcia Aguilar
import matplotlib.pyplot as plt #biblioteca grafica que nos va a permitir producir las graficas.
import math
import random #biblioteca que nos permite generar numeros al azar.
#trabajan en funciones
def validarADN(secuencia): #validaciones
    secuencia = secuencia.upper() #usamos solo mayusculas para evitar errores
    caracteresValidos = ["A", "C", "G", "T"] #definicion de los valores que vamos a permitir que se ingresen al programa
    for nucleotido in secuencia: #ciclo if que leera base por base e ira validando que esten en la lista de valores aceptables.
        if nucleotido not in caracteresValidos:
            return False
    return True


def analizarADN(secuencia): #conteos de cada una de las bases.
    secuencia = secuencia.upper()
    conteoA = secuencia.count('A') #contadores para cada una de las bases.
    conteoC = secuencia.count('C')
    conteoG = secuencia.count('G')
    conteoT = secuencia.count('T')

    complementoMap = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'} #hacemos una especie de diccionario con lo que fijamos cada letra con tu base correspondiente.
    cadenaComplementaria = ""

    for base in secuencia: #ya con el diccionario podemos generar el ciclo for para que recorra la cadena y llamamos al diccionario que generamos anteriormente.
        cadenaComplementaria += complementoMap[base]

    secuenciaARN = secuencia.replace('T', 'U') #generamos la cadena de arn, en donde reemplazamos cada tiamina por un uracilo.

    longitud = len(secuencia) #aqui guardamos la longitud de la secuencia ingresada.
    porcentajeGC = ((conteoG + conteoC) / longitud) * 100 if longitud > 0 else 0.0 #con esto calculamos el porcentajde de GC que hay dentro de la secuencia, con esto nos podemos dar idea de la estabilidad termica de la secuencia ya que la guanina y citosina formas un triple enlace a diferencia de adenina timina que solo forma un enlace doble.

    print("RESULTADOS DEL ANÁLISIS") #print de todo lo que calculamos o produjimos.
    print(f"Conteo = Adenina: {conteoA} | Citosina: {conteoC} | Guanina: {conteoG} | Timina: {conteoT}")
    print(f"Porcentaje Guanina-Citosina (GC): {porcentajeGC:.2f}%")
    print(f"Cadena Complementaria: {cadenaComplementaria}")
    print(f"ARN Mensajero Transcrito: {secuenciaARN}")

    return conteoA, conteoC, conteoG, conteoT, porcentajeGC #la devolucion de la informacion que se genero dentro de la funcion. Esta se puede usar posteriormente. Por cierto una funcion siempre te debe devolver un valor, ya sea numerico o booleano, pero siempre te regresa algo.


def traducirAminoacidos(secuenciaARN):
    codones = []
    cadenaProteina = ""
    
    return cadenaProteina


def simularMutacion(secuencia, probabilidad): #funcion con la cual simulamos las mutaciones dentro de una cadena.
    secuenciaMutada = "" #cadena vacia en donde se ira guardando la nueva cadena mutada.
    basesPosibles = ['A', 'C', 'G', 'T'] #de nuevo nuestra lista de posibles bases.
    
    for nucleotido in secuencia: #este for juega con la probabilidad que ingresa el usuario y los numeros producidos al azar para ir haciendole cambios a la cadena original.
        probGenerada = random.randint(1, 100)
        if probGenerada <= probabilidad:
            nuevaBase = random.choice(basesPosibles)
            secuenciaMutada += nuevaBase
        else:
            secuenciaMutada += nucleotido

    return secuenciaMutada #la funcion nos devuelve la cadena mutada.


def exportarYGraficar(secuencia, conteoA, conteoC, conteoG, conteoT, porcentajeGC): #aqui yo tampoco le se mucho pero es la biblioteca que genera los graficos.
    nombreArchivo = "Reporte Bioinformatico.txt" #nombre que tendra el archivo que se va a generar.
    try: #lo hacemos con un try para que no se caiga el programa en dado caso de que no se pueda generar el archivo.
        with open(nombreArchivo, "w") as archivo: #abre el archivo sobre el cual va a escribir y lo hace en modo de "escritura".
            archivo.write("REPORTE BIOINFORMATICO\n") #todo esto que dice "archivo.write" se escribira directamente sobre el archivo que se abrio, linea por lines, sera un documento .txt
            archivo.write(f"Secuencia Analizada: {secuencia}\n")
            archivo.write(f"Longitud Total: {len(secuencia)} pares de bases.\n\n")
            archivo.write(f"---CONTEO DE NUCLEOTIDOS---\n")
            archivo.write(f"Adenina (A): {conteoA}\n")
            archivo.write(f"Citosina (C): {conteoC}\n")
            archivo.write(f"Guanina (G): {conteoG}\n")
            archivo.write(f"Timina (T): {conteoT}\n")
            archivo.write(f"Estabilidad Termica (Contenido GC): {porcentajeGC:.2f}%\n")

        print(f"Éxito. Reporte textual guardado en '{nombreArchivo}'.") #maneja de manera segura el cierre del archivo.
    except Exception as e: #si no se puede abrir el archivo podremos ver un mensaje de error.
        print(f"Error. No se puede guardar el archivo: {e}")

    bases = ['Adenina (A)', 'Citosina (C)', 'Guanina (G)', 'Timina (T)'] 
    frecuencias = [conteoA, conteoC, conteoG, conteoT]
    colores = ["#E24AADEC", "#5072E3", "#B3F523", '#D0021B']

    plt.figure(figsize=(8, 5)) #tamaño del lienzo en pulgadas.
    plt.bar(bases, frecuencias, color=colores) #le decimos que genere una grafica de barras.
    plt.title('Frecuencia de Nucleotidos en la Secuencia', fontsize=14) #titulo del grafico
    plt.xlabel('Bases Nitrogenadas', fontsize=12)
    plt.ylabel('Cantidad', fontsize=12)
    print("Abriendo ventana de visualización gráfica...") #mansaje de aviso
    plt.show() #pedimos que se muestre el grafico en una ventana.


#modulo menu
def main():
    opcion = 0 #opciones del menu inicializadas en cero
    estadoDatosCorrectosCargados = False #es el estado que nos va a permitir corroborar que los datos ingresados sean correctos.
    secuenciaActual = "" #cadena vacia en donde se ira guardando la secuencia que ingrese el usuario

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
            opcion = int(input("Seleccione una opción: "))
        except ValueError:
            print("Error. Por favor ingrese un número.")
            continue

        if opcion == 1:
            entrada = input("Ingrese la secuencia de ADN: ")
            if validarADN(entrada):
                print("Secuencia válida.")
                secuenciaActual = entrada.upper()
                estadoDatosCorrectosCargados = True
            else:
                print("Error. La secuencia contiene caracteres no válidos.")
                estadoDatosCorrectosCargados = False

        elif opcion == 2:
            if estadoDatosCorrectosCargados:
                conteoA, conteoC, conteoG, conteoT, porcentajeGC = analizarADN(secuenciaActual)
            else:
                print("Se debe cargar una secuencia válida. Use la opción 1.")

        elif opcion == 3:
            if estadoDatosCorrectosCargados:
                probabilidad = int(input("Ingrese la probabilidad de mutación (1-100%): "))
                if 1 <= probabilidad <= 100:
                    secuenciaResultado = simularMutacion(secuenciaActual, probabilidad)
                    print(f"Secuencia Original: {secuenciaActual}")
                    print(f"Secuencia Mutada:   {secuenciaResultado}")
                else:
                    print("Error. Ingrese un valor entre 1 y 100.")
            else:
                print("Primero debe cargar una secuencia válida.")

        elif opcion == 4:
            if estadoDatosCorrectosCargados:
                conteoA, conteoC, conteoG, conteoT, porcentajeGC = analizarADN(secuenciaActual)
                exportarYGraficar(secuenciaActual, conteoA, conteoC, conteoG, conteoT, porcentajeGC)
            else:
                print("No hay datos para exportar.")

        elif opcion == 5:
            print("Saliendo del programa.")

        else:
            print("Opción no válida. Intente de nuevo.")


if __name__ == "__main__":
    main()

