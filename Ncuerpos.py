# -*- coding: utf-8 -*-
"""
Primera edición subido el 23 de abril de 2018

@author: eddybrandon
"""

import numpy as np
from scipy.constants import G

AU = (149597870700)     # Definiendo una UA
ao=1.2*10**(-10)         #Constante a_o

class Mond():
    nombre='ninguno'    #Será dado después.
    masa=0.0
    x=y=0.0
    vx=vy=0.0
    r=0.0
        
    #Aceleración a_ki de un cuerpo k sobre un cuerpo i:
    def aceleracion(self, otro):        #Aceleración k-i
        xki=self.x-otro.x
        yki=self.y-otro.y
        zki=self.z-otro.z
        rki=np.sqrt(xki**2+yki**2+zki**2) #Módulo radio r_i
        mk=otro.masa            #Masa del otro cuerpo.
        lmki=np.sqrt(G*mk/ao)   #Longitud característica.
        Xki=lmki/rki            #Cantidad adimensional.
        #funcion=Xki**2    #Función f(\chi) Newton (Utilizado para hacer los tests).
        funcion=Xki*(1+Xki+Xki**2+Xki**3)/(1+Xki+Xki**2)    #Función f(\chi) MOND
        ax=-ao*funcion*xki/rki
        ay=-ao*funcion*yki/rki
        az=-ao*funcion*zki/rki
        return ax, ay, az
        
        
def EulerN(astros):
    counter=0
    dt=24*3600   #Tamaño de paso = 1 día
    NI=5000      #Número de iteraciones
    arch={}
    for cuerpoi in astros:
        arch[cuerpoi.nombre]=open("Newton{}.dat".format(cuerpoi.nombre), "w")  #Archivos de texto de salida.
        
    while counter<NI:
        counter+=1
        accel={}    #Aceleración del cuerpo i
        for cuerpoi in astros:      #Sumatoria de aceleraciones a_ki sobre todos los k \neq i
            axi=ayi=azi=0
            if cuerpoi is 'Sol':    #En esta versión se fijó la posición del Sol en el origen.
                axi=0
                ayi=0
                azi=0
            else:
                for cuerpok in astros:
                    if cuerpok is cuerpoi:
                        continue
                    axki, ayki, azki=cuerpoi.aceleracion(cuerpok)
                    axi+=axki
                    ayi+=ayki
                    azi+=azki
            accel[cuerpoi]=(axi, ayi, azi)
            
        #MÉTODO DE EULER
        for cuerpoi in astros:
            axi, ayi, azi=accel[cuerpoi]
            cuerpoi.vx+=axi*dt      #Algoritmo v(n+1)=v(n)+a*dt
            cuerpoi.vy+=ayi*dt
            cuerpoi.vz+=azi*dt
            cuerpoi.x+=cuerpoi.vx*dt    # x(n+1)=x(n)+v(n)*dt
            cuerpoi.y+=cuerpoi.vy*dt
            cuerpoi.z+=cuerpoi.vz*dt
            escribir='{:>11.8f} {:>11.8f} {:>11.8f}'.format(cuerpoi.x/AU, cuerpoi.y/AU, cuerpoi.z/AU) #Datos x, y, z.
            arch[cuerpoi.nombre].write(escribir+"\n")     #Guarda las posiciones x_n, y_n, z_n en un archivo de texto para cada cuerpo.       

    for cuerpoi in astros:
        arch[cuerpoi.nombre].closed

#Lee las condiciones iniciales de cada planeta            
def condiciones(astros):
    for cuerpoi in astros:
        file=open('datos{}.dat'.format(cuerpoi.nombre), 'r')
        cuerpoi.masa = float(file.readline())       #Masa en [kg]
        cuerpoi.x = float(file.readline())*1000     #Posición en [km]
        cuerpoi.y = float(file.readline())*1000
        cuerpoi.z = float(file.readline())*1000
        cuerpoi.vx = float(file.readline())*1000    #Velocidad en [km/s]
        cuerpoi.vy = float(file.readline())*1000
        cuerpoi.vz = float(file.readline())*1000
        file.close        
        
def main():
    #Condiciones iniciales:
    
    sol = Mond()       
    sol.nombre = 'Sol'

    mercurio = Mond()
    mercurio.nombre = 'Mercurio'

    venus = Mond()
    venus.nombre = 'Venus'

    tierra = Mond()
    tierra.nombre = 'Tierra'

    marte = Mond()
    marte.nombre = 'Marte'
  
    jupiter = Mond()
    jupiter.nombre = 'Jupiter'
 
    saturno = Mond()
    saturno.nombre = 'Saturno'

    urano = Mond()
    urano.nombre = 'Urano'
    
    neptuno = Mond()
    neptuno.nombre = 'Neptuno'
    
    body = Mond()   #Cometa o Asteroide
    body.nombre = '16R2'
       
    astros=[sol, mercurio, venus, tierra, marte, jupiter, saturno, urano, neptuno, body]
    condiciones(astros)	        #Condiciones iniciales.
    EulerN(astros)              #Aplicación del método de Euler.
    
if __name__ == '__main__':
    main()
   
