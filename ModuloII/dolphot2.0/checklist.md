# PRACTICA MODULO II

### Culls de Rabdury-Smith 2011 [[10.1088/0067-0049/195/2/18](https://ui.adsabs.harvard.edu/link_gateway/2011ApJS..195...18R/doi:10.1088/0067-0049/195/2/18)]

* [X]  Cull3* (Sparse fields)
  * [X]  −0.06 < sharpnessF606W + sharpnessF814W < 1.30
  * [X]  crowdingF606W + crowdingF814W < 0.16
  * [X]  S/NF606W|F814W > 5.0
* [X]  Object type = 1
* [X]  error flag <= 2 en cada filtro

### Resultados

* [X]  Corregir por extincion galactica:  Monachesi+2016 $\rightarrow$ $A_V = 0.19$

  * [X]  From Schlafly & Finkbeiner, 2011 [10.1088/0004-637X/737/2/103](https://ui.adsabs.harvard.edu/link_gateway/2011ApJ...737..103S/doi:10.1088/0004-637X/737/2/103) - [NED database](https://irsa.ipac.caltech.edu/cgi-bin/bgTools/nph-bgExec)
  * [X]  f606w: 0.198 mags
  * [X]  f814w: 0.122 mags
    Se aplican los valores a cada valor de estrellas obtenido
* [X]  CMD

  Para fabricar el CMD se hace match en un radio de 80 pixeles (esto considerando la escala de pixel de la camara del hubble) para que las estrellas que esten en el diagrama correspondan en RA y DEC.
* [X]  Magnitud de la TRGB: 24.8, obtenido de la distribucion de las estrellas en el hexagono que corresponde al centro del CMD. La funcion de seaborn jointplot, entrega estos analisis estadisticos.

  Estimar $t_j$ :

  Partimos de la definicion entregada en clases de

  $t_j L_T B(t) = N_j$
  $t_j = \frac{N_j}{L_T B(t)}$
  donde reemplazamos con:  $B(t) = \frac{b(t)}{L_T(t)}$

  quedando $\rightarrow t_j = \frac{N_j}{b(t)}$ , con $b(t) = \phi(M_{TRGB})|\dot{M}_{TRGB}| $

  reemplazando:

  $t_j = \frac{N_{j}}{\phi(M_{TRGB}) |\dot{M}_{TRGB}|}  = \frac{N_{j}}{AM_{TRGB}^{-(1+x)} |\dot{M}_{TRGB}|}  $

  Dado que la perdida de masa y se estima usando la luminosidad, la cual podemos obtener desde los datos de las isocronas, finalmente $t_j$ quedaria:

  $t_j = \frac{N_{j}}{AM_{TRGB}^{-(1+x)} L^2Z}$

  con $N_j$ la cantidad de estrellas en el bin elegido (en este caso lo eligiriamos a partir de el jointplot) $x=1.5$ para utilizar una IMF de Salpeter y la luminosidad y Z extraida de los valores de la isocrona. Dado que el ajuste no fue satisfactorio, decidi no calcularlo.
* [X]  Ajustar isocrona para edad y metalicidad: Usando BasTi y los datos de [Monachesi+2013](https://iopscience.iop.org/article/10.1088/0004-637X/766/2/106), notamos que no es posible ajustar la isocrona a parte de lo que se observa en las estrellas, por lo cual, se plotea la isocrona mas correcta desde la literatura en naranjo, pero notamos que no se acerca a las estrellas de nuestro CMD. El principal problema esque la poblacion que se esta observando pareciera ser mucho mas metalica de lo que se puede acceder usando BasTi, y ademas, pareciera haber mas de una sola poblacion estelar. Aqui entramos a debatir acerca de la degeneracion edad metalicidad, y mas aun considerando que el filtro por el cual pasaron dejo pocas estrellas. No se descarta que hubieran errores asociados al calculo de la fotometria
* [X]  Nivel de completitud y poblacion del halo: Bajo, debido a los parametros utilizados. Pareciera ser que la poblacion estelar del halo, esta compuesta por mas de una sola poblacion, con diversas metalicidades, pero a las cuales no podemos profundizar debido a los limites instrumentales. Ademas, la region presenta un alto nivel de crowding, ya que luego de la correccion por el mismo que se propone en Radbury-Smit+2011, se pierden muchas fuentes detectadas.

### Imagenes adjuntas con los resultados:

* [X]  CMDs.jpg
* [X]  ischrone.jpg
* [X]  spatialdistribution.jpg
* [X]  trgbestimation.jpg
