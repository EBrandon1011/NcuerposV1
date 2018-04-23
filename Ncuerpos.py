# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 14:17:28 2018

@author: brandon
"""

import numpy as np
from scipy.constants import G

AU = (149597870700)     # Definiendo una ua
ao=1.2*10**(-10)         #Constante a_o

class Mond():
    nombre='ninguno'    #será dado después
    masa=0.0
    x=y=0.0
    vx=vy=0.0
    r=0.0
        
    def aceleracion(self, otro):        #Aceleración k-i
        xki=self.x-otro.x
        yki=self.y-otro.y
        rki=np.sqrt(xki**2+yki**2)
        #ri=np.sqrt(self.x**2+self.y**2)     #Módulo radio r_i
        mk=otro.masa
        lmki=np.sqrt(G*mk/ao)   #Longitud característica
        Xki=lmki/rki
        funcion=Xki**2    #Función f(\chi) Newton
        #funcion=Xki*(1+Xki+Xki**2+Xki**3)/(1+Xki+Xki**2)    #Función f(\chi) MOND
        ax=-ao*funcion*xki/rki
        ay=-ao*funcion*yki/rki
        #angi=np.arctan(self.y/self.x)     #Ángulo de r_i
        return ax, ay
        
        
def iteraciones(astros):
    counter=0
    dt=24*3600  #1 día
    arch={}
    for cuerpoi in astros:
        arch[cuerpoi.nombre]=open("Newton{}.dat".format(cuerpoi.nombre), "w")  #Archivos de texto
        
    while counter<400:
        counter+=1
        accel={}    #Aceleración sobre el cuerpo i
        for cuerpoi in astros:      #Sumatoria aceleraciones sobre k!=i
            axi=0
            ayi=0
            if cuerpoi is 'Sol':
                axi=0
                ayi=0
            else:
                for cuerpok in astros:
                    if cuerpok is cuerpoi:
                        continue
                    axki, ayki=cuerpoi.aceleracion(cuerpok)
                    axi+=axki
                    ayi+=ayki
            accel[cuerpoi]=(axi, ayi)
                
        for cuerpoi in astros:
            axi, ayi=accel[cuerpoi]
            cuerpoi.vx+=axi*dt      #Algoritmo v(n+1)=v(n)+a*dt
            cuerpoi.vy+=ayi*dt
            cuerpoi.x+=cuerpoi.vx*dt
            cuerpoi.y+=cuerpoi.vy*dt
            r=np.sqrt(cuerpoi.x**2+cuerpoi.y**2)
            ang=np.arctan(cuerpoi.y/cuerpoi.x)
            escribir='{:>11.8f} {:>11.8f}'.format(r/AU, ang) #Datos r, theta
            print(cuerpoi.nombre +': '+ escribir)
            arch[cuerpoi.nombre].write(escribir+"\n")            

    for cuerpoi in astros:
        arch[cuerpoi.nombre].closed
                        
def main():
    sol = Mond()        #Sol en el origen
    sol.nombre = 'Sol'
    sol.masa = 1.98892 * 10**30
    sol.x=0.0
    sol.vy=0.0

    tierra = Mond()
    tierra.nombre = 'Tierra'
    tierra.masa = 5.9742 * 10**24
    tierra.x = -1*AU
    tierra.vy = 29.783 * 1000  
    
    venus = Mond()
    venus.nombre = 'Venus'
    venus.masa = 4.8685 * 10**24
    venus.x = 0.723 * AU
    venus.vy = -35.02 * 1000
    
    iteraciones([sol, tierra, venus])
    
if __name__ == '__main__':
    main()
    
    
    
