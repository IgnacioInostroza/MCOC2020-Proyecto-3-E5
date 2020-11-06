from matplotlib.pylab import *
from matplotlib import cm


#definir geo

a=1.
b=1.
Nx=30
Ny=30

dx=b/Nx
dy=b/Ny

if dx!=dy:
    print("Error! dx!=dy")

h=dx

Puntos_a_analizar=[[],[],[]]

#funcion de conveniencia para calcular las coordenadas del punto (i,j)
coords = lambda i,j:(dx*i,dy*j)

x,y =coords(4,2)
print("x: ", x)
print("y: ", y)

def imshowbien(u):
    imshow(u.T[Nx::-1,:],cmap=cm.coolwarm, interpolation='bilinear')
    cbar=colorbar(extend='both',cmap=cm.coolwarm)
    ticks = arange(0,35,5)
    ticks_Text=["{}°".format(deg) for deg in ticks]
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(ticks_Text)
    clim(0,30)

    xlabel('a')
    ylabel('b')
    xTicks_N=arange(0,Nx+1,3)
    yTicks_N=arange(0,Ny+1,3)
    xTicks=[coords(i,0)[0] for i in xTicks_N]
    yTicks=[coords(0,i)[1] for i in yTicks_N]
    xTicks_Text=["{0:.2f}".format(tick) for tick in xTicks]
    yTicks_Text=["{0:.2f}".format(tick) for tick in yTicks]
    xticks(xTicks_N,xTicks_Text,rotation='vertical')
    yticks(yTicks_N,yTicks_Text)
    margins(0.2)
    subplots_adjust(bottom=0.15)
    #show()


u_k = zeros((Nx + 1, Ny + 1),dtype=double)
u_km1 = zeros((Nx + 1, Ny + 1),dtype=double)

#condiciones iniciales
u_k[:,:] = 20.  #20 grados inicial en todas partes

#parametros problema(hierro)
dt = 0.01 #s
K = 79.5 #m**2/s
c = 450.
rho = 7800.
alpha = K * dt / (c*rho * dx**2)

#Informar cosas interesantes
print(f"dt = {dt}")
print(f"dx = {dx}")
print(f"K = {K}")
print(f"c = {c}")
print(f"rho = {rho}")
print(f"alpha = {alpha}")

#Loop en el tiempo
minuto = 60.
hora = 3600.
dia = 24*3600

dt = 1*minuto
dnext_t = 0.5 * hora
next_t = 0
framenum = 0
T = 1*dia
Days = 1*T #cuantos dias quiero simular

def truncate(n,decimals=0):
    multiplier = 10**decimals
    return int(n*multiplier) / multiplier

#Loop en el tiempo
for k in range(int32(Days/dt)):
    t = dt*(k+1)
    dias = truncate(t/dia,0)
    horas = truncate((t-dias*dia)/hora,0)
    minutos = truncate((t-dias*dia-horas*hora)/minuto,0)
    titulo = "k = {0:05.0f}".format(k) + " t = {0:02.0f}d {1:02.0f}h {2:02.0f}m ".format(dias,horas,minutos)
    print(titulo)

    #CB esenciales, se repiten en cada iteracion
    u_k[0, :] = 20. #Borde izq
    u_k[-1, :] = 0.  # Borde der
    u_k[:, 0] = 20.  # Borde inferior
    u_k[:, -1] = 0.  # Borde sup

    #u_k[-1, :] = u_k[-2, :] - 10 * dx # Borde der, gradiente -10
    # (f(x+h)+f(x))/dx = algo

    #loop en el espacio i = 1...... n-1
    for i in range(1,Nx):
        for j in range(1,Ny):

            #algoritmo de diferencias finidas 2-D para difusion

            #laplaciano
            nabla_u_k = (u_k[i-1, j] + u_k[i+1, j] + u_k[i, j-1] + u_k[i, j+1] - 4 * u_k[i, j]) / h**2

            #forward Euler
            u_km1[i, j] = u_k[i, j]+alpha*nabla_u_k

    #avanzar la solucion a k+1
    u_k = u_km1

    #CB denuevo, para asegurar cumplimiento
    u_k[0, :] = 20. #Borde izq
    u_k[-1, :] = 0.  # Borde der
    u_k[:, 0] = 20.  # Borde inferior
    u_k[:, -1] = 0.  # Borde sup

    #u_k[-1, :] = u_k[-2, :] - 10 * dx # Borde der, gradiente -10
    # (f(x+h)+f(x))/dx = algo

    #Grafico en d_next
    if t>next_t:

        Puntos_a_analizar[0].append(u_k[int((Nx+1)/2),int((Ny+1)/2)])
        Puntos_a_analizar[1].append(u_k[int((Nx+1)/2),int(3*(Ny+1)/4)])
        Puntos_a_analizar[2].append(u_k[int(3*(Nx+1)/4),int(3*(Ny+1)/4)])

        figure(1)
        imshowbien(u_k)
        title(titulo)
        savefig("Ejemplo/frame_{0:04.0f}.png".format(framenum))
        framenum += 1
        next_t += dnext_t
        close(1)



x_x=np.linspace(dt,Days,len(Puntos_a_analizar[0]))


plt.plot(x_x,Puntos_a_analizar[0],label=" punto 1 (A/2, B/2 )")
plt.plot(x_x,Puntos_a_analizar[1],label=" punto 2 (A/2, 3B/4 )")
plt.plot(x_x,Puntos_a_analizar[2],label=" punto 3 (3A/4, 3B/4 )")
plt.legend(loc=1)
plt.title("Evolucion de temperatura en puntos")
plt.ylabel("Temeratura, $T$  [°C]")
plt.xlabel("Tiempo $t$ [segundos]")

plt.show()