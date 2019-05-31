Primero se carga inic_hmm.m
Ahí se ve que hay 6 emisones de cadenas de markov (X1...X6)

Cada uno de esos viene acompañado con un STi que marca que estado fue recorriendo la cadena cuando generó cada uno de estos cositos.

Las cadenas de markov están pensadas para que empiecen y nunca terminen (ergódicas) no nos sirve a nosotros porque las cosas están descriptas por cadenas de markov, hacen alguna cosa PERO EN ALGÚN MOMENTO van a un estado al cual no salen más.
Por lo tanto los estados ST tienen un estado inicial 1 y un estado final 5.
Correspondiente con el ST tiene que haber X con longitud igual a la longitud de ST-2 porque los estados 1 y 5 no emiten observaciones.

En habla tenemos un espacio de 39 dimensiones (39 estados).


Uno de los objetivos que vamos a tener es ver qué cadena de markov que está cargada (hmm_i) generó las muestras Xi.

hmm -> A (matriz de transición
-> b_j(x) = gaussiana

aij da la proba depasar de i en t-1 a j en t.

Ahora agarramos el X2 y vemos qué probabilidad tiene.

p(X,Z) = p(x1,z1) p(x2,z2) ... p(x9,z9)

Suponiendo que sólo depende del estado anterior (Sup de markov) =>
p(X,S) = p(s1=2)p(x1|s1=2) p(s2=2|s1=2)p(x2|s2=2) ... p(s9=4|s9=4)p(x9|s9=4)

Los p(x|s)=bj(x) son de los parámetros del modelo y p(s1=2|s)=aij sale de la primer fila de la matriz.


