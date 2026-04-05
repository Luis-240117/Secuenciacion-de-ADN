#Sofia Alejandra Beltran
#Luis Francisco Garcia Aguilar
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

#modulo del menu
def main():
    opcion=0
    estadoDatosCorrectosCargados=False
    secuenciaActual=""
    while opcion != 5:
        print("\n==========================")
        print("Sistema Life-Sequence")
        print("==========================")
        print("1. Cargar y Validar Secuencia de ADN")
        print("2. Analizar Secuencia Actual")
        print("3. Simular Mutaciones Aleatorias")
        print("4. Exportar Reporte y Graficar")
        print("5. Salir")