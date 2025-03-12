## **1. Introducción**

El estudio de la **dinámica molecular** es una herramienta clave en
bioinformática y más concretamente en modelado molecular. Esto permite
analizar el comportamiento y estabilidad de biomoléculas en condiciones
específicas. En este trabajo, se ha realizado una simulación de dinámica
molecular del tripéptido **ARV (Alanina-Arginina-Valina)** para estudiar
sus propiedades dinámicas y conformacionales.

Se han llevado a cabo simulaciones a **298 K** (condiciones
fisiológicas) y **400 K** (condiciones térmicas extremas) para evaluar
cómo la temperatura afecta la estabilidad estructural y dinámica del
tripéptido. Además, se incluye un análisis avanzado que abarca la
distribución de velocidades, gráficos de Ramachandran y la
identificación de conformaciones estructurales con mayor probabilidad.

## **2. Objetivos**

### **2.1. Preparar el sistema molecular:**

-   Visualizar el tripéptido ARV.

-   Definir la caja de simulación y solvatación.

### **2.2. Simulaciones de dinámica molecular:**

-   Realizar simulaciones en el ensamble **NVT** (Volumen, Número de
    partículas y Temperatura constante) a **298 K y 400 K**.

-   Registrar la evolución de energía, temperatura, geometría molecular
    y velocidades atómicas.

### **2.3. Análisis conformacional:**

-   Analizar propiedades geométricas como distancias de enlace, ángulos
    y dihedros.

-   Generar gráficos de Ramachandran para identificar las conformaciones
    más probables.

### **2.4.Comparación de resultados:**

-   Comparar los valores obtenidos en las simulaciones con datos
    teóricos disponibles en la literatura científica.

-   Evaluar cómo la temperatura afecta las propiedades dinámicas y
    estructurales del tripéptido.

### **2.5. Documentación y visualización:**

-   Elaborar gráficos claros y precisos utilizando `gnuplot` en `bash`.

-   Presentar los resultados en un informe en formato `.pdf`. Se realiza
    a partir de un RMarkdown que permite exportar a .pdf con LaTex y a
    otros formatos adecuados para ser publicado en GitHub.

## **3. Metodología y resultados**

### **3.1. Preparar el sistema molecular:**

-   Se crea un directorio de trabajo para resolver la tarea.

<!-- -->

    mkdir tarea

-   Se accede al directorio donde se encuentra el tripéptido y se copia
    el archivo necesario.

<!-- -->

    cd MM #nos situamos en el directorio donde está copiada for-students
    mkdir 1-build #creamos el directorio de trabajo para esta parte
    cd for-students/tripeptides #nos situamos en el directorio donde están copiados los tripéptidos
    ls #vemos el listado de tripéptidos
      aca  afa  aia  ama  aqa  ard  arg  ark  arn  ars  arw  ata  aya
      ada  aga  aka  ana  ara  are  arh  arl  arp  art  ary  ava
      aea  aha  ala  apa  arc  arf  ari  arm  arq  arv  asa  awa
    cp arv/arv.pdb /home/alumno08/MM/tarea/1-build #copiamos el directorio del tripéptido arv 
    ls -l #nos aseguramos de que se ha copiado en el directorio de trabajo

Ahora se va a usar el archivo PDB dentro del directorio arv para generar
los archivos necesarios para la simulación, creando los archivos de
topología necesarios para ello con el comando `pdb2gmx` desde el
directorio `1-build`.

    gmx pdb2gmx -f arv.pdb -o arv.gro -p arv.top -ter

**- Campo de fuerzas.**

GROMACS muestra una lista de campos de fuerzas disponibles y se elige el
8.

    8: CHARMM27 all-atom force field (CHARM22 plus CMAP for proteins)

**- Modelo de agua.**

Después de seleccionar el campo de fuerzas, se elige un modelo de agua.

    1: TIP3P   TIP 3-point, recommended

Una vez terminado, GROMACS genera los archivos `arv.gro` y `arv.top`.

    ls

    #Se obserevan los siguientes:
     arv.gro  arv.pdb  arv.top  posre.itp

Se visualiza `arv.gro`:

       62
        1ACE    CH3    1  -0.041   0.015   0.015
        1ACE   HH31    2  -0.053  -0.057   0.095
        1ACE   HH32    3  -0.076   0.111   0.049
        1ACE   HH33    4  -0.101  -0.017  -0.069
        1ACE      C    5   0.105   0.024  -0.024
        1ACE      O    6   0.166  -0.075  -0.066
        2ALA      N    7   0.166   0.137  -0.015
        2ALA     HN    8   0.118   0.219   0.019
        2ALA     CA    9   0.307   0.144  -0.054
        2ALA     HA   10   0.350   0.043  -0.050
        2ALA     CB   11   0.313   0.197  -0.198
        2ALA    HB1   12   0.258   0.132  -0.268
        2ALA    HB2   13   0.268   0.298  -0.207
    ...

Se visualiza `arv.top`.

    ;       This is a standalone topology file
    ;
    ;       Created by:
    ;                           :-) GROMACS - gmx pdb2gmx, 2016.4 (-:
    ;
    ;       Executable:   /usr/local/gromacs/bin/gmx
    ;       Data prefix:  /usr/local/gromacs
    ;       Working dir:  /home/alumno08/MM/tarea/1-build
    ;       Command line:
    ;         gmx pdb2gmx -f arv.pdb -o arv.gro -p arv.top -ter
    ;       Force field was read from the standard GROMACS share directory.
    ;

    ; Include forcefield parameters
    #include "charmm27.ff/forcefield.itp"

    [ moleculetype ]
    ; Name            nrexcl
    Protein             3
    ...

A partir de los archivos `arv.gro` y `arv.top`, se pueden interpretar
varios aspectos de la simulación molecular en GROMACS:

`arv.gro` es un archivo de coordenadas utilizado por GROMACS para
almacenar la posición de los átomos en el sistema. En él se observa que
hay 62 átomos en la estructura y la presencia de algunos residuos como
ACE (Acetil), que es un grupo bloqueante utilizado para evitar cargas
artificiales en los extremos de la cadena polipeptídica, y el residuo
ALA (Alanina), que indica que la simulación involucra una proteína.

`arv.top`define la topología del sistema, incluyendo fuerzas, tipos de
átomos, enlaces, restricciones, etc.

**- Caja de simulación.**

La caja de simulación define el espacio tridimensional en el que se
desarrollará la dinámica molecular. Se ejecuta el siguiente comando para
definir una caja cúbica de 3.0 nm por lado:

    Command line:
      gmx editconf -f arv.gro -o arv-box.gro -bt cubic -box 3.0 3.0 3.0

    Read 62 atoms
    Volume: 1.13989 nm^3, corresponds to roughly 500 electrons
    No velocities found
        system size :  1.359  0.665  1.261 (nm)
        diameter    :  1.525               (nm)
        center      :  0.640  0.289  0.213 (nm)
        box vectors :  1.359  0.665  1.261 (nm)
        box angles  :  90.00  90.00  90.00 (degrees)
        box volume  :   1.14               (nm^3)
        shift       :  0.860  1.211  1.287 (nm)
    new center      :  1.500  1.500  1.500 (nm)
    new box vectors :  3.000  3.000  3.000 (nm)
    new box angles  :  90.00  90.00  90.00 (degrees)
    new box volume  :  27.00               (nm^3)

`Expansión del volumen`

La caja inicial tenía un volumen de 1.14 nm³, lo que era insuficiente
para la simulación en condiciones realistas.

La nueva caja tiene un volumen de 27.00 nm³, permitiendo la adición de
solvente y evitando efectos de frontera (interacciones no físicas entre
imágenes periódicas de la molécula).

`Recentrado de la molécula`

La molécula estaba originalmente centrada en (0.640, 0.289, 0.213) nm.
Ahora está ubicada en (1.500, 1.500, 1.500) nm, lo que la posiciona en
el centro de la caja.

`Ángulos y vectores de la caja`

Se mantiene una caja ortorrómbica (90° en los tres ángulos). Los
vectores de la caja se ajustaron para reflejar el nuevo tamaño cúbico de
3.0 nm por lado.

`Falta de velocidades`

El mensaje *“No velocities found”* indica que este archivo
(`arv-box.gro`) sólo contiene posiciones atómicas, sin información de
velocidades iniciales. Esto es normal en este paso, ya que las
velocidades se asignarán más adelante en la simulación, normalmente a
través de un algoritmo de distribución de Maxwell-Boltzmann cuando se
inicie la dinámica molecular.

**- Solvatar la caja con agua.**

    gmx solvate -cp arv-box.gro -cs spc216.gro -o arv-solv.gro -p arv.top

En el archivo `arv-box-solv.gro` se muestran las coordenadas de las
moléculas de agua que se añadieron durante la etapa de solvatación.

Se visualiza `arv-box-solv.gro`.

    Great Red Oystrich Makes All Chemists Sane
     2636
        1ACE    CH3    1   0.819   1.226   1.302
        1ACE   HH31    2   0.807   1.154   1.382
        1ACE   HH32    3   0.784   1.322   1.336
        1ACE   HH33    4   0.759   1.194   1.218
        1ACE      C    5   0.965   1.235   1.263
        1ACE      O    6   1.026   1.136   1.221
        2ALA      N    7   1.026   1.348   1.272
        2ALA     HN    8   0.978   1.430   1.306
        2ALA     CA    9   1.167   1.355   1.233
        2ALA     HA   10   1.210   1.254   1.237
        2ALA     CB   11   1.173   1.408   1.089
        2ALA    HB1   12   1.118   1.343   1.019
        2ALA    HB2   13   1.128   1.509   1.080
        2ALA    HB3   14   1.276   1.415   1.052
        2ALA      C   15   1.245   1.441   1.330
        2ALA      O   16   1.203   1.550   1.370
        3ARG      N   17   1.360   1.399   1.372
    ...

-   Convertir el archivo solvatado a formato PDB.

<!-- -->

    gmx editconf -f arv-box-solv.gro -o arv-box-solv.pdb
    pymol arv-box-solv.pdb

<figure>
<img src="/home/alumno08/MM/imagen4.png"
alt="Representación del tripéptido ARV solvatado en Pymol" />
<figcaption aria-hidden="true">Representación del tripéptido ARV
solvatado en Pymol</figcaption>
</figure>

Se muestra la estructura del sistema tras la solvatación con agua, lo
que indica que se ha completado con éxito la etapa de añadir solvente a
la caja de simulación en GROMACS. El péptido ARV se visualiza en verde,
en el centro de la caja. Está rodeado por numerosas moléculas de agua
(oxígeno en rojo, hidrógeno en blanco), asegurando que la molécula se
encuentra en un entorno acuoso.

A continuación, se crea un directorio nuevo para este apartado y se
copian en él los ficheros requeridos.

    cd ..
    mkdir 2-equilibration
    cd 2-equilibration
    cp ../1-build/arv-box-solv.gro .
    cp ../1-build/arv.top .

### **3.2. Neutralización y equilibrado:**

**- A 298K.**

Se copia también el archivo `equiNVT.mdp` y se modifica cambaindo los
`nsteps` a 800000 para que al multiplicarse por los fs de 400 ps. En
este caso la temperatura se deja a 280K.

    nano equiNVT.mdp

    #Al ejecutar lo de arriba se obtiene:
    integrator = md ; leap-frog integrator
    dt = 0.0005 ; 0.5 fs
    nsteps = 800000 ; 400 ps #aquí se hace el cambio
    nstxout = 200 ; Save coordinates every 0.1 ps
    Tcoupl = v-rescale
    gen_vel = yes ; Should be the first equilibration
    gen_temp = 280.0 ; Temperature to generate corresponding Maxwell distribution
    gen_seed = 12345678 ; Random seed

    gmx grompp -f equiNVT.mdp -c arv-box-solv.gro -p arv.top -o arv-a.tpr

    #Al ejecutar lo de arriba se obtiene:
    NOTE 1 [file arv.top, line 613]:
      System has non-zero total charge: 1.000000

Como nuestra molécula tiene una carga total de +1 (*System has non-zero
total charge: 1.000000*) se le añade un Cl para neutralizarla.

    gmx genion -s arv-a.tpr -p arv.top -o arv-neutralized.gro -nn 1

    Select a group: 13
    Selected 13: 'SOL'
    Number of (3-atomic) solvent molecules: 858

    Processing topology
    Replacing 1 solute molecules in topology file (arv.top)  by 0 NA and 1 CL ions.

    Back Off! I just backed up arv.top to ./#arv.top.1#
    Using random seed -125359542.
    Replacing solvent molecule 706 (atom 2180) with CL

Este fragmento de salida muestra el ion añadido para neutralizar la
carga del sistema tras la solvatación en GROMACS. En muchas simulaciones
de dinámica molecular, se requiere que el sistema tenga carga neta cero
para evitar efectos artificiales en las interacciones electrostáticas.
En este caso, se selecciona el grupo SOL (agua), que en GROMACS
representa las moléculas de solvente y también se observa que el número
total de moléculas de agua en la caja es 858.

    gmx grompp -f equiNVT.mdp -c arv-neutralized.gro -p arv.top -o arv.tpr

    #Al ejecutar lo de arriba se obtiene:
    Analysing residue names:
    There are:     5    Protein residues
    There are:   857      Water residues
    There are:     1        Ion residues
    Analysing Protein...
    Analysing residues not classified as Protein/DNA/RNA/Water and splitting into groups...
    Number of degrees of freedom in T-Coupling group System is 7899.00
    Determining Verlet buffer for a tolerance of 0.005 kJ/mol/ps at 298 K
    Calculated rlist for 1x1 atom pair-list as 1.106 nm, buffer size 0.006 nm
    Set rlist, assuming 4x4 atom pair-list, to 1.100 nm, buffer size 0.000 nm

    gmx editconf -f arv-neutralized.gro -o arv.pdb
    pymol arv.pdb

<figure>
<img src="/home/alumno08/MM/imagen3.png"
alt="Representación del tripéptido ARV solvatado en y el ion Cl en Pymol" />
<figcaption aria-hidden="true">Representación del tripéptido ARV
solvatado en y el ion Cl en Pymol</figcaption>
</figure>

Se muestra el sistema después de la solvatación y adición de iones en la
caja de simulación. Como antes, se observa al péptido ARV representado
en verde con algunos átomos resaltados en azul y rodeado por moléculas
de agua (oxígeno en rojo, hidrógeno en blanco). El punto verde aislado
en la parte superior corresponde al ión Cl⁻ agregado en la etapa de
neutralización. En la simulación, este ion interactuará con la
biomolécula y el solvente, estabilizando la carga del sistema.

Para la ejecución, hay que copiarse el archivo `submit_eck.sh` al
directorio `2-equilibration`.

    cp /home/alumno08/MM/for-students/equilibration-mdp/submit_eck.sh .

`--chdir=$HOME/arv_simulation/2-equilibration`: el directorio donde está
la simulación. Hay que adaptarlo al `path` que interese.

    #!/bin/bash
    #
    #SBATCH -p eck-q
    #SBATCH --chdir=/home/alumno08/MM/tarea/2-equilibration
    #SBATCH -J equilibrado
    #SBATCH --cpus-per-task=1

    date
    gmx mdrun -deffnm arv -c arv.g96 -nt 1
    date

    sbatch submit_eck.sh

Se visualiza el archivo de salida `slurm-43311.out`.

    starting mdrun 'Protein in water'
    800000 steps,    400.0 ps.

    Writing final coordinates.

                   Core t (s)   Wall t (s)        (%)
           Time:     1512.091     1512.091      100.0
                     (ns/day)    (hour/ns)
    Performance:       22.856        1.050

El comando `gmx mdrun` ejecuta la simulación de producción tras la fase
de minimización de energía y equilibración.

En su `output`se observa que se ha corrido una dinámica molecular de 400
picosegundos (ps). El número total de pasos de integración fue 800,000,
lo que indica que el tamaño del paso de integración fue de 0.5
femtosegundos (fs). Tiempo total de ejecución: 1512.091 segundos (~25.2
minutos).

Se visualiza `arv.g96`.

    TITLE
    Protein in water
    END
    POSITION
        1 ACE   CH3        1    1.255372286    0.938619494    1.031705856
        1 ACE   HH31       2    1.196609974    0.911362350    0.939939380
        1 ACE   HH32       3    1.360876203    0.982189476    1.012639165
        1 ACE   HH33       4    1.270382881    0.840989769    1.092324972
        1 ACE   C          5    1.181761026    1.030441880    1.126240849
        1 ACE   O          6    1.161441326    1.142890692    1.093451619
        2 ALA   N          7    1.152526259    0.983615279    1.250377297
        2 ALA   HN         8    1.168643713    0.887859046    1.282865405
        2 ALA   CA         9    1.094711661    1.076678276    1.342497110
        2 ALA   HA        10    1.020556808    1.144909978    1.293710232
        2 ALA   CB        11    1.007464170    1.004780889    1.448202968
        2 ALA   HB1       12    1.065615535    0.928417027    1.502732158
        2 ALA   HB2       13    0.972558439    1.066741109    1.536105633
        2 ALA   HB3       14    0.920297503    0.954771221    1.399188280
        2 ALA   C         15    1.207374334    1.167288065    1.404257774
    ...

Con este archivo se obtienen los residuos a los que pertenecen los
átomos, así como el nombre átomo (CH3, HN, CB, etc.), su ID y las
coordenadas espaciales en nanómetros (nm).

Se visualiza `arv.log`.

    Initial temperature: 297.697 K

    Started mdrun on rank 0 Sat Mar  8 11:36:46 2025
               Step           Time
                  0        0.00000

       Energies (kJ/mol)
               Bond          Angle            U-B    Proper Dih.  Improper Dih.
        6.10355e+03    1.47996e+03    7.25045e+01    5.02421e+01    1.07150e-01
          CMAP Dih.          LJ-14     Coulomb-14        LJ (SR)  Disper. corr.
       -1.74058e+01    1.91861e+01   -1.07927e+03    3.59702e+04   -1.55883e+02
       Coulomb (SR)   Coul. recip.      Potential    Kinetic En.   Total Energy
       -3.82144e+04    1.04825e+03    5.27705e+03    1.19618e+04    1.72388e+04
      Conserved En.    Temperature Pres. DC (bar) Pressure (bar)
        1.72388e+04    3.64266e+02   -9.58703e+01    4.57637e+04

               Step           Time
                200        0.10000

       Energies (kJ/mol)
               Bond          Angle            U-B    Proper Dih.  Improper Dih.
        5.24368e+03    3.57990e+03    3.29951e+02    2.18057e+02    1.93595e+01
          CMAP Dih.          LJ-14     Coulomb-14        LJ (SR)  Disper. corr.
       -2.90784e+01    5.31018e+01   -9.74305e+02    7.79787e+03   -1.55883e+02
       Coulomb (SR)   Coul. recip.      Potential    Kinetic En.   Total Energy
       -4.18570e+04    3.67834e+02   -2.54065e+04    2.23378e+04   -3.06873e+03
      Conserved En.    Temperature Pres. DC (bar) Pressure (bar)
        1.58929e+04    6.80243e+02   -9.58703e+01   -7.07080e+03

    ...

Se muestra la evolución de energías, temperatura y presión a lo largo
del tiempo.

**- A 400K.**

Se repiten los pasos seguidos en el apartado de neutralización a 280K.
Para ello, se crea un directorio nuevo de trabajo y se copian en él los
archivos que se han usado para la neutralización a 280K.

Se modifica el archivo `equiNVT.mdp`, cambiando la temperatura se a
400K.

    ref_t = 400
    gen_temp = 400  

    gmx grompp -f equiNVT.mdp -c arv-box-solv.gro -p arv.top -o arv-a.tpr

    #Al ejecutar lo de arriba se obtiene:
    NOTE 1 [file arv.top, line 613]:
      System has non-zero total charge: 1.000000

Como nuestra molécula tiene una carga total de +1
(`System has non-zero total charge: 1.000000`) se le añade un Cl para
neutralizarla.

    gmx genion -s arv-a.tpr -p arv.top -o arv-neutralized.gro -nn 1

    Select a group: 13
    Selected 13: 'SOL'
    Number of (3-atomic) solvent molecules: 858

    Processing topology
    Replacing 1 solute molecules in topology file (arv.top)  by 0 NA and 1 CL ions.

    Back Off! I just backed up arv.top to ./#arv.top.1#
    Using random seed -118192417.
    Replacing solvent molecule 484 (atom 1514) with CL

    gmx grompp -f equiNVT.mdp -c arv-neutralized.gro -p arv.top -o arv.tpr

    #Al ejecutar lo de arriba se obtiene:
    Analysing residue names:
    There are:     5    Protein residues
    There are:   857      Water residues
    There are:     1        Ion residues
    Analysing Protein...
    Analysing residues not classified as Protein/DNA/RNA/Water and splitting into groups...
    Number of degrees of freedom in T-Coupling group System is 7899.00
    Determining Verlet buffer for a tolerance of 0.005 kJ/mol/ps at 400 K
    Calculated rlist for 1x1 atom pair-list as 1.108 nm, buffer size 0.008 nm
    Set rlist, assuming 4x4 atom pair-list, to 1.100 nm, buffer size 0.000 nm

Se definen los componentes del sistema, los grupos de acoplamiento de
temperatura, y los parámetros de la lista de vecinos para las
interacciones de pares atómicos. Se observa que el péptido ahora tiene
857 residuos de agua y 1 ion Cl⁻ agregado en la etapa de neutralización
del sistema. Con este `output` se concluye que el sistema está
correctamente preparado con proteína, solvente e iones. También se han
identificado los grados de libertad, lo que permite un control preciso
de la temperatura y se ha calculado un buffer de Verlet adecuado,
optimizando la simulación para precisión y rendimiento.

    #!/bin/bash
    #
    #SBATCH -p eck-q
    #SBATCH --chdir=/home/alumno08/MM/tarea/2-equilibration/400k
    #SBATCH -J equilibrado
    #SBATCH --cpus-per-task=1

    date
    gmx mdrun -deffnm arv -c arv.g96 -nt 1
    date

    sbatch submit_eck.sh

Se visualiza el archivo de salida `slurm-43314.out`.

    starting mdrun 'Protein in water'
    800000 steps,    400.0 ps.

    Writing final coordinates.

                   Core t (s)   Wall t (s)        (%)
           Time:     1553.856     1553.856      100.0
                     (ns/day)    (hour/ns)
    Performance:       22.241        1.079

El comando `gmx mdrun` ejecuta la simulación de producción tras la fase
de minimización de energía y equilibración.

En su `output`se observa que se ha corrido una dinámica molecular de 400
picosegundos (ps). El número total de pasos de integración fue 800,000,
lo que indica que el tamaño del paso de integración fue de 0.5
femtosegundos (fs). Tiempo total de ejecución: 1553.856 segundos (≈ 25.9
minutos).

Se visualiza `arv.g96`.

    TITLE
    Protein in water
    END
    POSITION
        1 ACE   CH3        1    2.586185932    1.792161345    0.572300732
        1 ACE   HH31       2    2.575107813    1.825858831    0.468126804
        1 ACE   HH32       3    2.553363323    1.873129725    0.631377399
        1 ACE   HH33       4    2.694044113    1.762442827    0.578654230
        1 ACE   C          5    2.511486292    1.670684457    0.609407544
        1 ACE   O          6    2.391749144    1.675015569    0.633933008
        2 ALA   N          7    2.589120865    1.559114337    0.604635477
        2 ALA   HN         8    2.681972980    1.560197234    0.582716286
        2 ALA   CA         9    2.547148228    1.423461199    0.636758208
        2 ALA   HA        10    2.640731812    1.382522225    0.670892119
        2 ALA   CB        11    2.442467928    1.409833193    0.755548954
        2 ALA   HB1       12    2.460976362    1.478632450    0.840936303
        2 ALA   HB2       13    2.338587046    1.424482584    0.730739594
    ...

Con este archivo se obtienen los residuos a los que pertenecen los
átomos, así como el nombre átomo (CH3, HN, CB, etc.), su ID y las
coordenadas espaciales en nanómetros (nm).

Se visualiza `arv.log`.

    Initial temperature: 399.594 K

    Started mdrun on rank 0 Sat Mar  8 12:37:02 2025
               Step           Time
                  0        0.00000

       Energies (kJ/mol)
               Bond          Angle            U-B    Proper Dih.  Improper Dih.
        6.10398e+03    1.47989e+03    7.25045e+01    5.02421e+01    1.07150e-01
          CMAP Dih.          LJ-14     Coulomb-14        LJ (SR)  Disper. corr.
       -1.74058e+01    1.91861e+01   -1.07927e+03    3.59581e+04   -1.55883e+02
       Coulomb (SR)   Coul. recip.      Potential    Kinetic En.   Total Energy
       -3.82385e+04    1.07672e+03    5.26966e+03    1.53914e+04    2.06610e+04
      Conserved En.    Temperature Pres. DC (bar) Pressure (bar)
        2.06610e+04    4.68707e+02   -9.58703e+01    4.71343e+04

               Step           Time
                200        0.10000

        Energies (kJ/mol)
               Bond          Angle            U-B    Proper Dih.  Improper Dih.
        5.94883e+03    4.27513e+03    3.76830e+02    2.09207e+02    3.62529e+01
          CMAP Dih.          LJ-14     Coulomb-14        LJ (SR)  Disper. corr.
       -2.75945e+01    5.49567e+01   -9.50585e+02    7.75606e+03   -1.55883e+02
       Coulomb (SR)   Coul. recip.      Potential    Kinetic En.   Total Energy
       -4.08482e+04    3.88598e+02   -2.29364e+04    2.50785e+04    2.14213e+03
      Conserved En.    Temperature Pres. DC (bar) Pressure (bar)
        1.93800e+04    7.63706e+02   -9.58703e+01   -6.34986e+03

               Step           Time
                400        0.20000
    ...

Se muestra la evolución de energías, temperatura y presión a lo largo
del tiempo.

### **3.3. Simulación de producción:**

**- A 298K.**

    cd ..
    mkdir 3-run
    cd 3-run
    cp ../2-equilibration/arv.top .
    cp ../for-students/g96-equilibrated/arv.g96 . 
    cp ../for-students/run-mdp/runNVT.mdp .

    gmx grompp -f runNVT.mdp -c arv.g96 -p arv.top -o arv.tpr

Archivos de salida `arv.trr` (coordenadas, velocidades, fuerzas,
trayectoria de la simulación)y `arv.tpr` (estructura, masas atómicas,
parámetros de simulación).

    #!/bin/bash
    #
    #SBATCH -p eck-q
    #SBATCH --chdir=/home/alumno08/MM/tarea/3-run/280k
    #SBATCH -J equilibrado
    #SBATCH --cpus-per-task=1  

    date
    gmx mdrun -deffnm arv -c arv.g96 -nt 1
    date

    sbatch submit_eck.sh
    squeue

    gmx traj -f arv.trr -s arv.tpr -oxt cartoon.pdb
    #Select: 0 (System)
    #Output file: cartoon.pdb

    pymol cartoon.pdb

<figure>
<img src="/home/alumno08/MM/tarea/3-run/280k/mi_peptido.png"
alt="Representación del tripéptido en Pymol tras la simulación a 298K" />
<figcaption aria-hidden="true">Representación del tripéptido en Pymol
tras la simulación a 298K</figcaption>
</figure>

Se concluye que la simulación ha sido enviada y ejecutada correctamente
en el clúster. El péptido se representa en verde con enlaces resaltados
en amarillo, las moléculas de agua (oxígeno en rojo, hidrógenos en
blanco) rodeando la proteína y el Cl⁻ en verde disperso en el solvente.
Las cargas electrostáticas están representadas como “+” y “-” indicando
interacciones iónicas.

    gmx traj -f arv.trr -s arv.tpr -oxt cartoon.pdb
    #Select: 1 (Protein)
    #Output file: cartoon.pdb

    pymol cartoon.pdb

Con este comando se extrae la estructura sin moléculas de agua ni iones.

<figure>
<img src="/home/alumno08/MM/tarea/3-run/280k/mi_peptido2.png"
alt="Representación en Pymol de la conformación final del péptido tras la simulación a 298K" />
<figcaption aria-hidden="true">Representación en Pymol de la
conformación final del péptido tras la simulación a 298K</figcaption>
</figure>

**- A 400K.**

    cd ..
    mkdir 3-run
    cd 3-run
    cp ../2-equilibration/arv.top .
    cp ../for-students/g96-equilibrated/arv.g96 . 
    cp ../for-students/run-mdp/runNVT.mdp .

    gmx grompp -f runNVT.mdp -c arv.g96 -p arv.top -o arv.tpr

    #!/bin/bash
    #
    #SBATCH -p eck-q
    #SBATCH --chdir=/home/alumno08/MM/tarea/3-run/400k
    #SBATCH -J equilibrado
    #SBATCH --cpus-per-task=1  

    date
    gmx mdrun -deffnm arv -c arv.g96 -nt 1
    date

    sbatch submit_eck.sh
    squeue

    gmx traj -f arv.trr -s arv.tpr -oxt cartoon.pdb
    #Select: 0 (System)
    #Output file: cartoon.pdb

    pymol cartoon.pdb

    gmx traj -f arv.trr -s arv.tpr -oxt cartoon.pdb
    #Select: 1 (Protein)
    #Output file: cartoon.pdb

    pymol cartoon.pdb

### **3.4. Análisis:**

    cd ..
    mkdir 4-analysis
    cd 4-analysis
    cp ../3-run/arv.tpr .
    cp ../3-run/arv.g96 .
    cp ../3-run/arv.trr .  # Coordenadas
    cp ../3-run/arv.edr .  # Energías

**- A 298K.**

    gmx gyrate -f arv.trr -s arv.tpr -xvg none
    #Select 1 (Protein)

    gnuplot> plot 'gyrate.xvg'

<figure>
<img src="/home/alumno08/MM/tarea/4-analysis/280k/1.png"
alt="Representación de la evolución de giro del péptido a 298K" />
<figcaption aria-hidden="true">Representación de la evolución de giro
del péptido a 298K</figcaption>
</figure>

La gráfica representa la evolución del radio de giro (`Rg`) del péptido
durante la simulación de dinámica molecular. El radio de giro es una
medida del tamaño global de la proteína y refleja cambios en su
compactación o expansión a lo largo del tiempo. En ella se muestran
muestran fluctuaciones del Rg, lo que indica que el péptido está en
constante movimiento, alternando entre estados más compactos y más
extendidos.

El eje X representa el tiempo de simulación (en nanosegundos o
picosegundos).

El eje Y muestra el radio de giro en nanómetros (nm).

Al final de este apartado se interpreta y compara esta gráfica a los
289K con la correspondiente a 400K.

**Distancias de enlace.**

Distancia entre el primer carbono alfa (CA) y el primer carbono beta
(CB) del péptido.

**Identificación de los átomos relevantes.**

    gmx dump -s arv.tpr | grep "CA" #Para ver qué residuos son los primeros CA y CB
    gmx dump -s arv.tpr | grep "CB"

    #Esto nos devuelve:
    atom[8]={name="CA"} atom[18]={name="CA"} atom[42]={name="CA"}
    atom[10]={name="CB"} atom[20]={name="CB"} atom[44]={name="CB"}


    nano distances.ndx

    #Dentro del archivo:

    [CA_CB]
    8 10

**Cálculo de distancias a lo largo del tiempo.**

    gmx distance -f arv.trr -s arv.tpr -n distances.ndx -oall bond_distances.xvg -xvg none
    #Select group: 0
    #Ctrl+d

    #Se obtiene:
    Number of samples:  2001
      Average distance:   0.28016  nm
      Standard deviation: 0.01154  nm

`Distancia media entre CA y CB`: 0.280 nm

`Desviación estándar`: 0.0115 nm, lo que indica una variabilidad
relativamente baja en la distancia.

Esto sugiere que la distancia CA-CB se mantiene estable con ligeras
fluctuaciones a lo largo de la simulación.

**Visualizar los resultados.**

    less ca_cb_distance.xvg

    #Se observa:

          0.000    0.297
          0.001    0.294
          0.002    0.290
          0.003    0.285
          0.004    0.280
          0.005    0.277
          0.006    0.276
          0.007    0.277
          0.008    0.279
          0.009    0.282

Con este archivo se observa que la distancia comienza en 0.297 nm, pero
se estabiliza en torno a 0.280 nm. También se aprecia una ligera
reducción inicial, lo que podría indicar un ajuste estructural en la
conformación del péptido.

Posteriormente, las fluctuaciones son mínimas, lo que sugiere que la
distancia CA-CB está relativamente restringida y estable.

    gnuplot

    gnuplot> plot "ca_cb_distance.xvg" with linespoints

<figure>
<img src="/home/alumno08/MM/tarea/4-analysis/280k/2.png"
alt="Representación de la distancia de los CA y CB seleccionados a 298K" />
<figcaption aria-hidden="true">Representación de la distancia de los CA
y CB seleccionados a 298K</figcaption>
</figure>

Esta gráfica muestra la evolución de la distancia entre el carbono alfa
(CA) y el carbono beta (CB) del péptido a lo largo del tiempo, con datos
extraídos del archivo `ca_cb_distance.xvg`. En el eje X se representa el
tiempo (en ns o ps), mientras que el eje Y indica la distancia en
nanómetros (nm). La distancia CA-CB se mantiene estable con ligeras
fluctuaciones, reflejando la flexibilidad estructural del péptido,
aunque se observa una leve tendencia a la disminución, lo que podría
sugerir un pequeño reacomodo en la conformación inicial. No se presentan
cambios conformacionales abruptos, lo que indica que la estructura del
péptido se mantiene relativamente rígida en este intervalo de
simulación. Al final de este apartado, se interpreta y compara esta
gráfica a 298K con la correspondiente a 400K.

**Siguiente paso.**

GROMACS tiene una herramienta para encontrar el contacto más cercano
entre átomos.

Este análisis tiene como objetivo encontrar el oxígeno de carbonilo (CO)
más cercano a los carbonos alfa (CA) y beta (CB) del péptido.

    gmx dump -s arv.tpr | grep "atom"

    #Se obtiene los siguientes residuos:

    4 (C), 5 (O)
    14 (C), 15 (O)
    38 (C), 39 (O)
    54 (C), 55 (O)

Elegimos el CO de los residuos 4 y 5.

El siguiente comando mide la evolución de la distancia entre los átomos
especificados en `distances.ndx` a lo largo de la simulación.

    #Modificamos distances.ndx

    [CA_CB]
    8 10
    [CO]
    4 5

    gmx distance -f arv.trr -s arv.tpr -n distances.ndx -oall distances_CA_CB_CO.xvg -xvg none

    #Se obtiene:

    Analyzed 2001 frames, last time 2.000
    CO:
      Number of samples:  2001
      Average distance:   0.21377  nm
      Standard deviation: 0.00616  nm

La distancia promedio CO es 0.214 nm (2.14 Å).

Esta distancia sugiere que el grupo CO está cercano al CA y CB
seleccionados, lo cual es esperable en la estructura peptídica.

La desviación estándar es 0.00616 nm. Una baja variabilidad sugiere que
el contacto entre el CO y los carbonos seleccionados es estable durante
la simulación.

La proximidad del CO al CA y CB podría indicar una posible interacción
estabilizadora, por ejemplo, puentes de hidrógeno entre el grupo
carbonilo y los átomos cercanos o estructuras secundarias donde el CO
forma parte de un motivo estructural más grande (hojas β, hélices α).

<figure>
<img src="/home/alumno08/MM/tarea/4-analysis/280k/4.png"
alt="Representación de la distancia del CA, CB y CO seleccioandos a 298K" />
<figcaption aria-hidden="true">Representación de la distancia del CA, CB
y CO seleccioandos a 298K</figcaption>
</figure>

Esta gráfica muestra la evolución de la distancia entre los átomos **CA,
CB y CO** a lo largo del tiempo en la simulación de dinámica molecular,
con el eje **X** representando el tiempo (en ns o ps) y el eje **Y** la
distancia en nanómetros (nm). Inicialmente, la distancia ronda los 0.30
nm y presenta fluctuaciones constantes, observándose una leve
disminución hasta estabilizarse en torno a 0.27-0.28 nm en la segunda
mitad de la simulación, lo que indica que los átomos mantienen una
proximidad relativamente constante con cierta flexibilidad. Las
oscilaciones en la distancia (~0.02-0.03 nm) son esperadas en dinámica
molecular debido a los movimientos naturales de la estructura, y la
ausencia de cambios abruptos sugiere que no hay una reconfiguración
estructural significativa. La estabilidad de la distancia indica una
posible interacción persistente entre estos átomos, como un puente de
hidrógeno o interacciones estabilizadoras dentro de la estructura
secundaria, hipótesis que se analiza en detalle más adelante mediante el
cálculo de puentes de hidrógeno con `gmx hbond`.

**Ángulos de enlace.**

Se define un grupo que contenga CA, CB y CO del mismo residuo o residuos
cercanos.

    nano distances.ndx
    #Se actualiza el contenido

    [CA-CB-CO]
    8 10 4

    gmx angle -f arv.trr -n distances.ndx -ov angle_CA_CB_CO.xvg -xvg none


    #Se obtiene:
    Found points in the range from 16 to 49 (max 180)
     < angle >  = 29.9125
    < angle^2 > = 942.31
    Std. Dev.   = 6.89585

    gnuplot> plot "angle_CA_CB_CO.xvg" with linespoints

<figure>
<img src="/home/alumno08/MM/tarea/4-analysis/280k/3.png"
alt="Representación del ángulo formado por el grupo molecular CA-CB-CO a 298K" />
<figcaption aria-hidden="true">Representación del ángulo formado por el
grupo molecular CA-CB-CO a 298K</figcaption>
</figure>

Esta gráfica representa la evolución del ángulo formado entre los átomos
**CA, CB y CO** a lo largo del tiempo en la simulación de dinámica
molecular. En el eje **X** se encuentra el tiempo (en ns o ps), mientras
que en el eje **Y** se representa el valor del ángulo en grados. Se
observa que el ángulo oscila entre **15° y 50°**, con variaciones
significativas a lo largo de la trayectoria. Estas fluctuaciones
sugieren que la orientación relativa de los átomos experimenta cambios
constantes, lo cual podría indicar una flexibilidad estructural en la
región donde se encuentran. No se observa una tendencia clara a la
estabilización en un valor fijo, lo que sugiere que el ángulo sigue
evolucionando dinámicamente sin alcanzar una conformación estable a lo
largo del tiempo.

**Ángulos dihédricos.**

Este análisis se centra en la extracción y procesamiento de los
**ángulos Phi y Psi**.

    gmx rama -f arv.trr -s arv.tpr -xvg none -o rama_arv.xvg

    less rama_arv.xvg

    #Se obtiene:
    -79.9678  151.902  ALA-2
    -70.5712  143.844  ARG-3
    -132.615  -34.3638  VAL-4
    -80.3427  151.36  ALA-2
    -69.1537  143.378  ARG-3
    -132.95  -34.1079  VAL-4
    -80.6473  150.769  ALA-2
    -67.7175  142.933  ARG-3
    -133.234  -33.8205  VAL-4
    -80.8407  150.124  ALA-2
    -66.347  142.539  ARG-3

El comando genera un archivo `rama_arv.xvg` con los valores de los
ángulos Phi y Psi para cada residuo de interés en la simulación.

Cada línea del `output` representa los valores de Phi y Psi para los
residuos ALA-2, ARG-3 y VAL-4 en distintos frames de la simulación.

    awk '{print $1}' angle_CA_CB_CO.xvg | cat -n > phi_CA_CB_CO.dat
    awk '{print $2}' angle_CA_CB_CO.xvg | cat -n > psi_CA_CB_CO.dat

    wc -l phi_CA_CB_CO.dat
    2001 phi_CA_CB_CO.dat

    wc -l psi_CA_CB_CO.dat
    2001 psi_CA_CB_CO.dat

    join phi_CA_CB_CO.dat psi_CA_CB_CO.dat > angles_CA_CB_CO.dat

Finalmente, se organizan los valores individuales con `awk`, y se
combinan en un solo archivo con `join` para facilitar su análisis. Estos
datos pueden utilizarse para evaluar la estabilidad conformacional y la
flexibilidad de los residuos en la simulación.

<figure>
<img src="/home/alumno08/MM/tarea/4-analysis/280k/5.png"
alt="Representación del ángulo Phi a 298K" />
<figcaption aria-hidden="true">Representación del ángulo Phi a
298K</figcaption>
</figure>

<figure>
<img src="/home/alumno08/MM/tarea/4-analysis/280k/6.png"
alt="Representación del ángulo Psi a 298K" />
<figcaption aria-hidden="true">Representación del ángulo Psi a
298K</figcaption>
</figure>

La gráfica muestra la evolución del ángulo **Psi** entre los átomos
**CA, CB y CO** a lo largo del tiempo en la simulación de dinámica
molecular. En el eje **X** se representa el número de frames de la
simulación, mientras que en el eje **Y** se observa el valor del ángulo
en grados, con variaciones entre aproximadamente **15° y 50°**. Como en
los parámetros anteriores, se evidencia un comportamiento altamente
dinámico con oscilaciones continuas, lo que sugiere una flexibilidad
estructural considerable en la región donde se encuentran estos átomos.
No se observa una estabilización clara en un valor fijo, lo que indica
que el ángulo sigue fluctuando activamente a lo largo de la simulación,
posiblemente debido a cambios conformacionales locales.

<figure>
<img src="/home/alumno08/MM/tarea/4-analysis/280k/7.png"
alt="Representación del ángulo CA-CB-CO a 298K" />
<figcaption aria-hidden="true">Representación del ángulo CA-CB-CO a
298K</figcaption>
</figure>

La gráfica muestra la evolución del ángulo **CA-CB-CO** a lo largo del
tiempo en la simulación de dinámica molecular y presenta un
comportamiento muy similar al del otro ángulo analizado. Se observan
oscilaciones entre **15° y 50°**, sin una tendencia clara a la
estabilización, lo que indica una conformación flexible y dinámica en
esta región de la estructura. Al igual que en la otra gráfica, los
cambios en el ángulo sugieren una reorganización constante,
probablemente debido a interacciones transitorias o a la movilidad
natural del péptido.

**Enlaces de hidrógeno (H-bonds).**

Este análisis evalúa la presencia y estabilidad de puentes de hidrógeno
(H-bonds).

    gmx hbond -f arv.trr -s arv.tpr -a 120 -r 0.18 -noda -num num_hbond.xvg -xvg none -nonitacc -hx hx_hbond.xvg

    #Select 1 and 1 (Protein) for intramolecular
    #Select 1 and 13 (SOL) for intermolecular

    less num_hbond.xvg
    #Se obtiene:
             0           5           0
         0.001           5           0
         0.002           4           0
         0.003           4           0
         0.004           5           0
         0.005           4           0
         0.006           4           0
         0.007           3           0
         0.008           3           0
         0.009           3           0
    ...
    less hx_hbond.xvg
    #Se obtiene:
             0      0      0      0      0      0      0      0
         0.001      0      0      0      0      0      0      0
         0.002      0      0      0      0      0      0      0
         0.003      0      0      0      0      0      0      0
         0.004      0      0      0      0      0      0      0
         0.005      0      0      0      0      0      0      0
         0.006      0      0      0      0      0      0      0
    ...

Los resultados en `num_hbond.xvg` indican la cantidad de H-bonds en cada
frame, mientras que visualizando `hx_hbond.xvg`se deduce que no hay
hidrógenos que se mantengan involucrados de forma constante en puentes
de hidrógeno, lo que sugiere que los H-bonds detectados son transitorios
y no estables.

La proteína mantiene entre 3 y 5 puentes de hidrógeno internos durante
la simulación, pero no se detectan interacciones fuertes con el solvente
bajo los criterios definidos. Además, los H-bonds parecen transitorios,
ya que ningún hidrógeno mantiene una participación estable en estas
interacciones. Si se quisiera evaluar con más detalle la estabilidad de
estos puentes de hidrógeno, se podría analizar su persistencia a lo
largo del tiempo o modificar los criterios de radio y ángulo para
detectar interacciones más débiles.

**Velocidad y temperatura.**

    gmx traj -f arv.trr -s arv.tpr -xvg none -ot temp_arv.xvg

    #Select a group: 0

    gmx traj -f arv.trr -s arv.tpr -xvg none -ov -len
    #Select a group: 1

    less veloc.xvg

    #Se obtiene:
     0      -0.223168       -0.211114       0.00986779      0.30736 2.27404 0.292277        -0.500722       2.34679 2.01421 3.1197  -2.10053        4.26636 -0.0447923      -1.14785        2.54943 2.79627 0.132697        -0.669125       0.00546973      0.682178        0.479224        -0.225444       -0.290474       0.604033        -0.178657       1.33826 -0.353151       1.39556 0.773205        0.414245        -0.835374       1.21132 0.26341 0.688991        -0.570203       0.932322        -1.051  -2.07189        -1.30009        2.66225 0.894955        0.426097        -0.505292       1.11257 -0.218003       4.04986 -0.193359       4.06033 

El primer comando permite extraer la evolución de la temperatura a lo
largo del tiempo. Este análisis permite verificar si la temperatura se
mantiene estable o si hay fluctuaciones significativas que puedan
afectar la dinámica del sistema.

El segundo comando extrae información sobre las velocidades de los
átomos en la simulación. Cada línea de su `output` representa un frame
de la simulación seguido por los valores de velocidad en las tres
direcciones (x, y, z) para cada átomo. Estos valores permiten evaluar la
dinámica molecular del sistema, identificando regiones con alta o baja
movilidad atómica.

<figure>
<img src="/home/alumno08/MM/tarea/4-analysis/280k/8.png"
alt="Representación de la evolución de temperatura a 298K" />
<figcaption aria-hidden="true">Representación de la evolución de
temperatura a 298K</figcaption>
</figure>

La gráfica muestra la evolución de la temperatura del sistema a lo largo
del tiempo en la simulación de dinámica molecular, con los datos
extraídos del archivo `temp_arv.xvg`. En el eje **X** se representa el
tiempo en nanosegundos (ns), mientras que en el eje **Y** se observa la
temperatura en Kelvin (K), con fluctuaciones entre aproximadamente **280
K y 315 K**. Se aprecia una variabilidad considerable en la temperatura,
pero sin desviaciones extremas, lo que indica que el sistema mantiene un
control térmico razonable dentro de los límites esperados. Sin embargo,
las oscilaciones pueden estar relacionadas con la elección del
termostato en la simulación o con eventos estructurales en la proteína
que influyen en la distribución de la energía. Para un análisis más
detallado, se podría calcular la temperatura media y su desviación
estándar para verificar si se mantiene cerca del valor objetivo.

    # Cargar el archivo de temperatura
    temp_data <- read.table("/home/alumno08/MM/tarea/4-analysis/280k/temp_arv.xvg", header=FALSE)

    # Renombrar las columnas (asumiendo que la primera es el tiempo y la segunda es la temperatura)
    colnames(temp_data) <- c("Time", "Temperature")

    # Calcular la temperatura media y desviación estándar
    temp_mean <- mean(temp_data$Temperature)
    temp_sd <- sd(temp_data$Temperature)

    # Mostrar los resultados
    cat("Temperatura media:", temp_mean, "K\n")

    ## Temperatura media: 297.2643 K

    cat("Desviación estándar:", temp_sd, "K\n")

    ## Desviación estándar: 4.922913 K

<figure>
<img src="/home/alumno08/MM/tarea/4-analysis/280k/10.png"
alt="Representación de la evolución de velocidad a 298K" />
<figcaption aria-hidden="true">Representación de la evolución de
velocidad a 298K</figcaption>
</figure>

La gráfica muestra la evolución de la velocidad de un átomo específico a
lo largo del tiempo en la simulación de dinámica molecular, con datos
extraídos del archivo `veloc.xvg`. En el eje **X** se representa el
tiempo en nanosegundos (ns), mientras que en el eje **Y** se muestra la
magnitud de la velocidad en nanómetros por picosegundo (nm/ps). Se
observa una distribución dispersa de los valores de velocidad, con
oscilaciones en el rango de **0 a 1.8 nm/ps**, lo que sugiere una
movilidad considerable. La ausencia de una tendencia clara indica que la
velocidad varía de manera aleatoria, lo cual es característico de
sistemas en equilibrio térmico. Para una mejor interpretación, se podría
calcular la **media y la desviación estándar** de las velocidades, o
analizar su distribución con un histograma para verificar si sigue una
distribución normal esperada en un sistema térmico equilibrado.

    gnuplot> n=100 #number of intervals
    gnuplot> max=3. #max value
    gnuplot> min=-3. #min value
    gnuplot> width=(max-min)/n #interval width
    gnuplot> hist(x,width)=width*floor(x/width)+width/2.0
    gnuplot> set boxwidth width*0.9
    gnuplot> set style fill solid 0.5 # fill style
    gnuplot> plot 'veloc.xvg' u (hist($2,width)):(1.0) smooth freq w boxes lc rgb'red' notitle
    gnuplot> plot 'veloc.xvg' u (hist($3,width)):(1.0) smooth freq w boxes lc rgb'red' notitle

<figure>
<img src="/home/alumno08/MM/tarea/4-analysis/280k/11.png"
alt="Histograma de velocidades con hist($2,width) a 298K" />
<figcaption aria-hidden="true">Histograma de velocidades con
<code>hist($2,width)</code> a 298K</figcaption>
</figure>

<figure>
<img src="/home/alumno08/MM/tarea/4-analysis/280k/12.png"
alt="Histograma de velocidades con hist($3,width) a 298K" />
<figcaption aria-hidden="true">Histograma de velocidades con
<code>hist($3,width)</code> a 298K</figcaption>
</figure>

Ambos histogramas representan la distribución de las velocidades
atómicas en dos componentes diferentes de la simulación de dinámica
molecular, mostrando una forma de campana característica de una
**distribución normal**, lo que indica un comportamiento térmicamente
equilibrado. Sin embargo, existen diferencias notables: en la primera
gráfica, la distribución está más concentrada alrededor de **0 nm/ps**,
con un pico más pronunciado, lo que sugiere que la mayoría de los átomos
tienen velocidades moderadas en esa dirección; mientras que en la
segunda gráfica, el pico es más bajo y la dispersión es ligeramente
mayor, lo que indica una mayor variabilidad en esa componente de
velocidad. Ambas distribuciones cubren un rango similar de valores,
aproximadamente entre **-1.5 y 1.5 nm/ps**, pero la ligera diferencia en
su dispersión sugiere que una de las direcciones tiene más fluctuaciones
dinámicas que la otra, lo que podría estar relacionado con la
anisotropía del sistema o la interacción con el entorno.

    gmx energy -f arv.edr -s arv.tpr -xvg none

    -------------------------------------------------------------------
      1  Bond             2  Angle            3  U-B              4  Proper-Dih.
      5  Improper-Dih.    6  CMAP-Dih.        7  LJ-14            8  Coulomb-14
      9  LJ-(SR)         10  Disper.-corr.   11  Coulomb-(SR)    12  Coul.-recip.
     13  Potential       14  Kinetic-En.     15  Total-Energy    16  Conserved-En.
     17  Temperature     18  Pres.-DC        19  Pressure        20  Vir-XX
     21  Vir-XY          22  Vir-XZ          23  Vir-YX          24  Vir-YY
     25  Vir-YZ          26  Vir-ZX          27  Vir-ZY          28  Vir-ZZ
     29  Pres-XX         30  Pres-XY         31  Pres-XZ         32  Pres-YX
     33  Pres-YY         34  Pres-YZ         35  Pres-ZX         36  Pres-ZY
     37  Pres-ZZ         38  #Surf*SurfTen   39  T-System        40  Lamb-System

    15

    Last energy frame read 2000 time    2.000

    Statistics over 4001 steps [ 0.0000 through 2.0000 ps ], 1 data sets
    All statistics are over 2001 points (frames)

    Energy                      Average   Err.Est.       RMSD  Tot-Drift
    -------------------------------------------------------------------------------
    Total Energy               -26124.4         85    217.018   -397.378  (kJ/mol)

<figure>
<img src="/home/alumno08/MM/tarea/4-analysis/280k/13.png"
alt="Representación de la energía del sistema a 298K" />
<figcaption aria-hidden="true">Representación de la energía del sistema
a 298K</figcaption>
</figure>

La gráfica representa la evolución de la **energía total** del sistema a
lo largo del tiempo en la simulación de dinámica molecular, con datos
extraídos del archivo `energy.xvg`. En el eje **X** se encuentra el
tiempo en **nanosegundos (ns)**, mientras que en el eje **Y** se muestra
la energía total en **kJ/mol**. Se observa que la energía oscila entre
aproximadamente **9200 y 10300 kJ/mol**, con fluctuaciones constantes
que reflejan el intercambio de energía entre los diferentes grados de
libertad del sistema. Aunque la variabilidad es notable, el sistema no
muestra una tendencia clara de incremento o disminución significativa,
lo que sugiere que se ha alcanzado un estado de equilibrio dinámico.
Según los valores obtenidos en **GROMACS**, la **energía total
promedio** es de **-26124.4 kJ/mol**, con una desviación **RMSD de
217.018 kJ/mol** y una **deriva total de -397.378 kJ/mol**, lo que
indica una ligera pérdida de energía a lo largo de la simulación.

    -------------------------------------------------------------------
      1  Bond             2  Angle            3  U-B              4  Proper-Dih.
      5  Improper-Dih.    6  CMAP-Dih.        7  LJ-14            8  Coulomb-14
      9  LJ-(SR)         10  Disper.-corr.   11  Coulomb-(SR)    12  Coul.-recip.
     13  Potential       14  Kinetic-En.     15  Total-Energy    16  Conserved-En.
     17  Temperature     18  Pres.-DC        19  Pressure        20  Vir-XX
     21  Vir-XY          22  Vir-XZ          23  Vir-YX          24  Vir-YY
     25  Vir-YZ          26  Vir-ZX          27  Vir-ZY          28  Vir-ZZ
     29  Pres-XX         30  Pres-XY         31  Pres-XZ         32  Pres-YX
     33  Pres-YY         34  Pres-YZ         35  Pres-ZX         36  Pres-ZY
     37  Pres-ZZ         38  #Surf*SurfTen   39  T-System        40  Lamb-System

    14


    Back Off! I just backed up energy.xvg to ./#energy.xvg.1#
    Last energy frame read 2000 time    2.000

    Statistics over 4001 steps [ 0.0000 through 2.0000 ps ], 1 data sets
    All statistics are over 2001 points (frames)

    Energy                      Average   Err.Est.       RMSD  Tot-Drift
    -------------------------------------------------------------------------------
    Kinetic En.                 9765.34         26     160.16   -165.267  (kJ/mol)

El análisis de la **energía cinética** en la simulación de dinámica
molecular muestra un valor promedio de **9765.34 kJ/mol**, con una
desviación cuadrática media (**RMSD**) de **160.16 kJ/mol** y un error
estimado de **26 kJ/mol**, lo que indica fluctuaciones moderadas pero
esperadas en un sistema en equilibrio térmico. Además, la deriva total
de **-165.267 kJ/mol** sugiere una ligera disminución en la energía
cinética a lo largo del tiempo, posiblemente debido a la disipación de
energía controlada por el termostato utilizado en la simulación. Estas
variaciones son normales y reflejan la redistribución de energía entre
los grados de libertad del sistema.

**- A 400K.**

Se repite el anális anterior para la simulación a 400K.

    gmx gyrate -f arv.trr -s arv.tpr -xvg none
    #Select 1 (Protein)

    gnuplot> plot 'gyrate.xvg'

<figure>
<img src="/home/alumno08/MM/tarea/4-analysis/400k/1.png"
alt="Representación de la evolución de giro del péptido a 400K" />
<figcaption aria-hidden="true">Representación de la evolución de giro
del péptido a 400K</figcaption>
</figure>

La gráfica representa la evolución del **radio de giro (`Rg`)** del
péptido durante la simulación a **400K**, reflejando cambios en su
compactación y expansión a lo largo del tiempo. Comparando con la
simulación a **298K**, se observa que en ambos casos el **Rg fluctúa**,
lo que indica que el péptido mantiene su dinamismo estructural en ambas
temperaturas. Sin embargo, a **400K** las variaciones parecen más
pronunciadas y frecuentes, lo que sugiere una mayor flexibilidad
estructural debido al aumento de la temperatura. En este caso, la
proteína experimenta **transiciones más marcadas entre estados compactos
y extendidos**, lo que es consistente con un mayor movimiento molecular
a temperaturas elevadas. Esta diferencia podría estar relacionada con
una menor estabilidad estructural o con un posible inicio de
**desnaturalización térmica** en el péptido, dependiendo de la magnitud
de las fluctuaciones observadas.

**Distancias de enlace.**

Distancia entre el primer carbono alfa (CA) y el primer carbono beta
(CB) del péptido.

**Identificación de los átomos relevantes.**

    gmx dump -s arv.tpr | grep "CA" #Para ver qué residuos son los primeros CA y CB
    gmx dump -s arv.tpr | grep "CB"

    #Esto nos devuelve:
    atom[8]={name="CA"} atom[18]={name="CA"} atom[42]={name="CA"}
    atom[10]={name="CB"} atom[20]={name="CB"} atom[44]={name="CB"}


    nano distances.ndx

    #Dentro del archivo:

    [CA_CB]
    8 10

**Cálculo de distancias a lo largo del tiempo.**

    gmx distance -f arv.trr -s arv.tpr -n distances.ndx -oall bond_distances.xvg -xvg none
    #Select group: 0
    #Ctrl+d

    #Se obtiene:
      Number of samples:  1829
      Average distance:   0.21656  nm
      Standard deviation: 0.01290  nm
      
    less bond_distances.xvg

    #Se obtiene:
          0.000    0.203
          0.001    0.201
          0.002    0.200
          0.003    0.200
          0.004    0.201
          0.005    0.203
          0.006    0.204
          0.007    0.207
          ...

<figure>
<img src="/home/alumno08/MM/tarea/4-analysis/400k/2.png"
alt="Representación de la distancia de los CA y CB seleccionados a 400K" />
<figcaption aria-hidden="true">Representación de la distancia de los CA
y CB seleccionados a 400K</figcaption>
</figure>

La gráfica muestra la evolución de la **distancia entre el carbono alfa
(CA) y el carbono beta (CB)** del péptido a **400K**, con el eje **X**
representando el tiempo (en ns) y el eje **Y** la distancia en
**nanómetros (nm)**. En comparación con la simulación a **298K**, donde
la distancia se mantenía estable con ligeras fluctuaciones en torno a
**0.27-0.28 nm**, a **400K la distancia promedio es de 0.21656 nm**, con
una **desviación estándar de 0.01290 nm**, lo que indica una ligera
reducción en la separación entre estos átomos. Se observan
**fluctuaciones más pronunciadas**, lo que sugiere un aumento en la
flexibilidad estructural a mayor temperatura. Sin embargo, al igual que
en **298K**, no se presentan cambios conformacionales abruptos, lo que
indica que la estructura del péptido sigue siendo relativamente estable
en este intervalo de simulación, aunque nuevamente con mayor movilidad
molecular.

**Siguiente paso.**

GROMACS tiene una herramienta para encontrar el contacto más cercano
entre átomos.

Se define un grupo que contenga CA, CB y CO del mismo residuo o residuos
cercanos.

    gmx dump -s arv.tpr | grep "atom"

    #Se obtiene los siguientes residuos:

    4 (C), 5 (O)
    14 (C), 15 (O)
    38 (C), 39 (O)
    54 (C), 55 (O)

Elegimos el CO de los residuos 4 y 5.

    #Modificamos distances.ndx

    [CA_CB]
    8 10
    [CO]
    4 5

    gmx distance -f arv.trr -s arv.tpr -n distances.ndx -oall distances_CA_CB_CO.xvg -xvg none

    #Se obtiene:

    Analyzed 1829 frames, last time 1.829
    CO:
      Number of samples:  1829
      Average distance:   0.21362  nm
      Standard deviation: 0.00680  nm

<figure>
<img src="/home/alumno08/MM/tarea/4-analysis/400k/3.png"
alt="Representación de la distancia de los CA, CB y CO seleccionados a 400K" />
<figcaption aria-hidden="true">Representación de la distancia de los CA,
CB y CO seleccionados a 400K</figcaption>
</figure>

La gráfica muestra la evolución de la distancia entre los átomos **CA,
CB y CO** a lo largo del tiempo en la simulación de dinámica molecular a
**400K**. En comparación con la simulación a **298K**, donde la
distancia oscilaba alrededor de **0.27-0.28 nm**, en **400K la distancia
media se reduce a aproximadamente 0.216 nm**, con fluctuaciones más
marcadas a lo largo del tiempo. Aunque no se observan cambios abruptos
que indiquen una reconfiguración estructural drástica, la reducción de
la distancia y el incremento en la amplitud de las fluctuaciones podrían
reflejar **ajustes conformacionales y una posible modificación en la
estabilidad de las interacciones moleculares** dentro de la estructura.
Para un análisis más detallado, se evalua la persistencia de puentes de
hidrógeno mediante `gmx hbond`, ya que la menor distancia podría estar
relacionada con la formación o ruptura de interacciones estabilizadoras.

**Ángulos de enlace.**

Se define un grupo que contenga CA, CB y CO del mismo residuo o residuos
cercanos.

    nano distances.ndx
    #Se actualiza el contenido

    [CA-CB-CO]
    8 10 4

    gmx angle -f arv.trr -n distances.ndx -ov angle_CA_CB_CO.xvg -xvg none


    #Se obtiene:
    Found points in the range from 2 to 39 (max 180)
     < angle >  = 20.4177
    < angle^2 > = 487.157
    Std. Dev.   = 8.38296

    gnuplot> plot "angle_CA_CB_CO.xvg" with linespoints

<figure>
<img src="/home/alumno08/MM/tarea/4-analysis/400k/4.png"
alt="Representación del ángulo CA-CB-CO a 400K" />
<figcaption aria-hidden="true">Representación del ángulo CA-CB-CO a
400K</figcaption>
</figure>

La gráfica muestra la evolución del **ángulo CA-CB-CO** a lo largo del
tiempo en la simulación a **400K**, presentando un comportamiento
similar al observado a **298K**, pero con diferencias clave en la
magnitud y la frecuencia de las oscilaciones. Mientras que a **298K** el
ángulo fluctuaba entre **15° y 50°**, en esta simulación las
oscilaciones son más amplias y frecuentes, alcanzando valores más bajos
(alrededor de **5°**) y evidenciando una mayor variabilidad. Esto
sugiere que el aumento de la temperatura a **400K** ha incrementado la
**flexibilidad conformacional** de esta región del péptido, lo que
podría estar asociado a una mayor libertad de movimiento y
reorganización estructural. La falta de una tendencia clara a la
estabilización indica que el ángulo sigue fluctuando de manera dinámica,
posiblemente debido a interacciones transitorias o cambios en la
conformación del péptido inducidos por la temperatura.

**Ángulos dihédricos.**

    gmx rama -f arv.trr -s arv.tpr -xvg none -o rama_arv.xvg

    less rama_arv.xvg

    #Se obtiene:
    99.7554  158.644  ALA-2
    -147.872  141.913  ARG-3
    -75.2336  131.001  VAL-4
    98.9408  158.538  ALA-2
    -148.179  141.712  ARG-3
    -75.6868  129.837  VAL-4
    98.0417  158.44  ALA-2
    -148.319  141.591  ARG-3
    -76.4718  128.921  VAL-4
    97.0973  158.316  ALA-2
    -148.305  141.544  ARG-3
    ...

    awk '{print $1}' angle_CA_CB_CO.xvg | cat -n > phi_CA_CB_CO.dat
    awk '{print $2}' angle_CA_CB_CO.xvg | cat -n > psi_CA_CB_CO.dat

    wc -l phi_CA_CB_CO.dat
    1829 phi_CA_CB_CO.dat

    wc -l psi_CA_CB_CO.dat
    1829 psi_CA_CB_CO.dat

    join phi_CA_CB_CO.dat psi_CA_CB_CO.dat > angles_CA_CB_CO.dat

<figure>
<img src="/home/alumno08/MM/tarea/4-analysis/400k/5.png"
alt="Representación del ángulo Phi a 400K" />
<figcaption aria-hidden="true">Representación del ángulo Phi a
400K</figcaption>
</figure>

<figure>
<img src="/home/alumno08/MM/tarea/4-analysis/400k/6.png"
alt="Representación del ángulo Psi a 400K" />
<figcaption aria-hidden="true">Representación del ángulo Psi a
400K</figcaption>
</figure>

Mientras que en **298K** las variaciones del ángulo oscilaban entre
**15° y 50°**, en esta simulación a **400K** se observa que las
oscilaciones son más pronunciadas, alcanzando valores más bajos
(alrededor de **5°**), lo que indica un aumento en la flexibilidad
estructural. Al igual que en la simulación a **298K**, el ángulo no
presenta una tendencia clara hacia la estabilización, lo que confirma
que sigue fluctuando dinámicamente debido a la mayor energía térmica
presente en el sistema a **400K**.

<figure>
<img src="/home/alumno08/MM/tarea/4-analysis/400k/7.png"
alt="Representación del ángulo CA-CB-CO a 400K" />
<figcaption aria-hidden="true">Representación del ángulo CA-CB-CO a
400K</figcaption>
</figure>

La gráfica del **ángulo CA-CB-CO** a **400K** muestra oscilaciones más
amplias y frecuentes en comparación con **298K**, indicando una mayor
**flexibilidad estructural** debido al aumento de temperatura. Aunque el
ángulo sigue fluctuando sin estabilizarse, la mayor variabilidad sugiere
un incremento en la **movilidad conformacional**.

**Enlaces de hidrógeno (H-bonds).**

Este análisis evalúa la presencia y estabilidad de puentes de hidrógeno
(H-bonds).

    gmx hbond -f arv.trr -s arv.tpr -a 120 -r 0.18 -noda -num num_hbond.xvg -xvg none -nonitacc -hx hx_hbond.xvg

    #Select 1 and 1 (Protein) for intramolecular
    #Select 1 and 13 (SOL) for intermolecular

    less num_hbond.xvg

    #Se obtiene:
             0           0           0
         0.001           0           0
         0.002           0           0
         0.003           0           0
         0.004           0           0
         0.005           0           0
         0.006           0           0
         0.007           0           0
         0.008           0           0
         0.009           0           0
          0.01           0           0
         0.011           0           0
         0.012           1           0
        ...

    less hx_hbond.xvg
    #Se obtiene:
         0.057      0      0      0      0      0      0      0
         0.058      0      0      0      0      0      0      0
         0.059      0      0      0      0      0      0      0
          0.06      0      0      0      0      0      0      0
         0.061      0      0      0      0      0      0      0
         0.062      0      0      0      0      0      0      0
         0.063      0      0      0      0      0      0      0
         0.064      0      0      0      0      0      0      0
         0.065      0      0      0      0      0      0      0
         0.066      0      0      0      0      0      0      0
         0.067      0      0      0      0      0      0      0
         0.068      0      0      0      0      0      0      0
         0.069      0      0      0      0      0      0      0
        ...

A **400K**, los puentes de hidrógeno intramoleculares son escasos y
transitorios, lo que sugiere una disrupción estructural en comparación
con simulaciones a temperaturas más bajas como 298K, donde había una
mayor presencia de interacciones estabilizadoras.

    gmx traj -f arv.trr -s arv.tpr -xvg none -ot temp_arv.xvg

    #Select a group: 0

    gmx traj -f arv.trr -s arv.tpr -xvg none -ov -len
    #Select a group: 1

    less veloc.xvg

    #Se obtiene:
      0      -0.553119       -0.119066       -0.607265       0.829993        0.12516 1.09539 -2.77298        2.98412 1.46621 0.58409 2.41471 2.88475 1.03624 1.7385  1.30308 2.40711 0.0624985       0.0936952       0.317636        0.337013        -0.315851       0.326307        0.127426        0.471673        -0.192531       0.549627        0.459647        0.741912        0.771657        -3.3487 -0.213139       3.44307 -0.526927       0.20526 -1.30288        1.42031 -1.48797        -0.140763       -1.836  2.36744 1.17616 0.0721901       0.16804 1.19029 -1.20551        1.45106 -2.20746        2.90374 -2.78303        -0.0384135      -0.269188       2.79628 0.935413        -2.646  1.46685 3.1667  -0.762276       0.555381        0.500643        1.06778 -0.375983       -0.109886       -0.719887       0.819558        -0.75507        0.633363        0.708175        1.21359 -0.35363
      ...

Un **aumento de temperatura** favorece la movilidad de los átomos y
puede llevar a una pérdida de interacciones estabilizadoras, como se
observa en la disminución de los puentes de hidrógeno.

La mayor energía cinética a 400K aumenta la **movilidad atómica**, lo
que puede facilitar cambios conformacionales y una disrupción de la
estructura secundaria.

<figure>
<img src="/home/alumno08/MM/tarea/4-analysis/400k/11.png"
alt="Representación de la evolución de la temperatura del sistema a 400K" />
<figcaption aria-hidden="true">Representación de la evolución de la
temperatura del sistema a 400K</figcaption>
</figure>

La gráfica muestra la evolución de la **temperatura del sistema** a
**400K** a lo largo del tiempo en la simulación de dinámica molecular,
con datos extraídos del archivo `temp_arv.xvg`. En el eje **X** se
representa el tiempo en **nanosegundos (ns)** y en el eje **Y** la
**temperatura en Kelvin (K)**, con fluctuaciones entre aproximadamente
**380 K y 415 K**. Comparado con la simulación a **298K**, donde la
temperatura oscilaba entre **280 K y 315 K**, se observa una mayor
variabilidad térmica, lo que sugiere que el sistema está sometido a una
mayor agitación molecular. A pesar de estas fluctuaciones, la
temperatura se mantiene dentro del rango esperado, indicando que el
**control térmico** funciona correctamente. Para una evaluación más
precisa, sería útil calcular la **temperatura media y su desviación
estándar**, permitiendo verificar si se mantiene estable en torno a
**400K**.

    # Cargar el archivo de temperatura
    temp_data <- read.table("/home/alumno08/MM/tarea/4-analysis/400k/temp_arv.xvg", header=FALSE)

    # Renombrar las columnas (asumiendo que la primera es el tiempo y la segunda es la temperatura)
    colnames(temp_data) <- c("Time", "Temperature")

    # Calcular la temperatura media y desviación estándar
    temp_mean <- mean(temp_data$Temperature)
    temp_sd <- sd(temp_data$Temperature)

    # Mostrar los resultados
    cat("Temperatura media:", temp_mean, "K\n")

    ## Temperatura media: 397.8735 K

    cat("Desviación estándar:", temp_sd, "K\n")

    ## Desviación estándar: 5.632015 K

<figure>
<img src="/home/alumno08/MM/tarea/4-analysis/400k/8.png"
alt="Representación de la evolución de la velocidad del sistema a 400K" />
<figcaption aria-hidden="true">Representación de la evolución de la
velocidad del sistema a 400K</figcaption>
</figure>

La gráfica muestra la evolución de la **velocidad de un átomo
específico** en la simulación a **400K**, con datos extraídos del
archivo `veloc.xvg`. En comparación con la simulación a **298K**, donde
las velocidades fluctuaban en un rango de **0 a 1.8 nm/ps**, en **400K**
se aprecian valores ligeramente más dispersos, con algunos puntos
alcanzando velocidades superiores a **2.0 nm/ps**, lo que indica una
mayor movilidad atómica debido al incremento de la temperatura.

    gnuplot> n=100 #number of intervals
    gnuplot> max=3. #max value
    gnuplot> min=-3. #min value
    gnuplot> width=(max-min)/n #interval width
    gnuplot> hist(x,width)=width*floor(x/width)+width/2.0
    gnuplot> set boxwidth width*0.9
    gnuplot> set style fill solid 0.5 # fill style
    gnuplot> plot 'veloc.xvg' u (hist($2,width)):(1.0) smooth freq w boxes lc rgb'red' notitle
    gnuplot> plot 'veloc.xvg' u (hist($3,width)):(1.0) smooth freq w boxes lc rgb'red' notitle

<figure>
<img src="/home/alumno08/MM/tarea/4-analysis/400k/9.png"
alt="Histograma de velocidades con hist($2,width) a 400K" />
<figcaption aria-hidden="true">Histograma de velocidades con
<code>hist($2,width)</code> a 400K</figcaption>
</figure>

<figure>
<img src="/home/alumno08/MM/tarea/4-analysis/400k/10.png"
alt="Histograma de velocidades con hist($3,width) a 400K" />
<figcaption aria-hidden="true">Histograma de velocidades con
<code>hist($3,width)</code> a 400K</figcaption>
</figure>

Las gráficas muestran la distribución de las **velocidades atómicas** en
la simulación a **400K**, con histogramas obtenidos a partir del archivo
`veloc.xvg`. En ambas, el eje **X** representa la velocidad (en nm/ps),
mientras que el eje **Y** muestra la frecuencia de ocurrencia de cada
valor de velocidad.

Al comparar entre **298K y 400K** se deduce que, en ambas temperaturas,
las velocidades siguen una **distribución normal centrada en 0 nm/ps**,
característica de un sistema en equilibrio térmico. Sin embargo, a
**400K** se observa una mayor dispersión en los valores de velocidad, lo
que indica un incremento en la movilidad atómica respecto a **298K**.
Además, la menor altura del pico en **400K** sugiere una distribución
más extendida de las velocidades, reflejando el aumento de la **energía
cinética** debido a la mayor temperatura.

    gmx energy -f arv.edr -s arv.tpr -xvg none

    -------------------------------------------------------------------
      1  Bond             2  Angle            3  U-B              4  Proper-Dih.
      5  Improper-Dih.    6  CMAP-Dih.        7  LJ-14            8  Coulomb-14
      9  LJ-(SR)         10  Disper.-corr.   11  Coulomb-(SR)    12  Coul.-recip.
     13  Potential       14  Kinetic-En.     15  Total-Energy    16  Conserved-En.
     17  Temperature     18  Pres.-DC        19  Pressure        20  Vir-XX
     21  Vir-XY          22  Vir-XZ          23  Vir-YX          24  Vir-YY
     25  Vir-YZ          26  Vir-ZX          27  Vir-ZY          28  Vir-ZZ
     29  Pres-XX         30  Pres-XY         31  Pres-XZ         32  Pres-YX
     33  Pres-YY         34  Pres-YZ         35  Pres-ZX         36  Pres-ZY
     37  Pres-ZZ         38  #Surf*SurfTen   39  T-System        40  Lamb-System

    15

    Last energy frame read 2000 time    2.000

    Statistics over 4001 steps [ 0.0000 through 2.0000 ps ], 1 data sets
    All statistics are over 2001 points (frames)

    Energy                      Average   Err.Est.       RMSD  Tot-Drift
    -------------------------------------------------------------------------------
    Total Energy               -17425.8        8.3    322.409    3.58843  (kJ/mol)

<figure>
<img src="/home/alumno08/MM/tarea/4-analysis/400k/13.png"
alt="Representación de la energía del sistema a 400K" />
<figcaption aria-hidden="true">Representación de la energía del sistema
a 400K</figcaption>
</figure>

En la segunda gráfica se observa un comportamiento similar a **298K**,
pero con una **mayor amplitud en las fluctuaciones energéticas**,
reflejando una mayor movilidad molecular y mayor agitación térmica
debido al aumento de temperatura. La energía total oscila en un rango
más amplio (**12200 a 14200 kJ/mol**), lo que es consistente con el
incremento de la energía cinética del sistema. A pesar de la
variabilidad, el sistema parece mantenerse en equilibrio, sin cambios
drásticos en la tendencia de la energía total. Sin embargo, la mayor
dispersión de los valores podría indicar una mayor flexibilidad
estructural o un posible inicio de desestabilización en la conformación
del péptido a esta temperatura.

    -------------------------------------------------------------------
      1  Bond             2  Angle            3  U-B              4  Proper-Dih.
      5  Improper-Dih.    6  CMAP-Dih.        7  LJ-14            8  Coulomb-14
      9  LJ-(SR)         10  Disper.-corr.   11  Coulomb-(SR)    12  Coul.-recip.
     13  Potential       14  Kinetic-En.     15  Total-Energy    16  Conserved-En.
     17  Temperature     18  Pres.-DC        19  Pressure        20  Vir-XX
     21  Vir-XY          22  Vir-XZ          23  Vir-YX          24  Vir-YY
     25  Vir-YZ          26  Vir-ZX          27  Vir-ZY          28  Vir-ZZ
     29  Pres-XX         30  Pres-XY         31  Pres-XZ         32  Pres-YX
     33  Pres-YY         34  Pres-YZ         35  Pres-ZX         36  Pres-ZY
     37  Pres-ZZ         38  #Surf*SurfTen   39  T-System        40  Lamb-System

    14


    Back Off! I just backed up energy.xvg to ./#energy.xvg.1#
    Last energy frame read 2000 time    2.000

    Statistics over 4001 steps [ 0.0000 through 2.0000 ps ], 1 data sets
    All statistics are over 2001 points (frames)

    Energy                      Average   Err.Est.       RMSD  Tot-Drift
    -------------------------------------------------------------------------------
    Kinetic En.                 13134.8        2.5    207.518    12.7479  (kJ/mol)

El análisis de la **energía cinética** en la simulación indica un valor
promedio de **13134.8 kJ/mol**, con una desviación **RMSD de 207.518
kJ/mol** y un error estimado de **2.5 kJ/mol**, lo que sugiere
fluctuaciones moderadas pero esperadas en un sistema en equilibrio
térmico. La **deriva total de 12.7479 kJ/mol** es mínima, indicando que
la energía cinética se mantiene estable a lo largo de la simulación.
Estos valores reflejan el impacto del aumento de temperatura en la
movilidad molecular, con una mayor energía cinética en comparación con
simulaciones a temperaturas más bajas.

### **3.5. Tarea extra 1:**

    cd ..
    mkdir 3-run
    cd 3-run
    cp ../2-equilibration/arv.top .
    cp ../280k/g96-equilibrated/arv.g96 . 
    cp ../for-students/run-mdp/runNVT.mdp .

    nano runNVT.mdp

    #Se modifican estos parámetros
    nsteps = 1000000 ; 500 ps
    nstxout = 20 ; Save coordinates every 10 fs
    nstvout = 20 ; Save velocities every 10 fs
    nstlog = 20 ; Update log every 10 fs
    nstenergy = 20 ; Save energies every 10 fs

    sbatch submit_eck.sh

Se visualiza el archivo de salida `slurm-43314.out`.

    starting mdrun 'Protein in water'
    1000000 steps,    500.0 ps.

    Writing final coordinates.

    Back Off! I just backed up arv.g96 to ./#arv.g96.1#

                   Core t (s)   Wall t (s)        (%)
           Time:     1914.780     1914.780      100.0
                             31:54
                     (ns/day)    (hour/ns)
    Performance:       22.561        1.064

Se visualiza `arv.g96`.

    TITLE
    Protein in water
    END
    POSITION
        1 ACE   CH3        1    1.122357965    1.024833083    2.110969782
        1 ACE   HH31       2    1.187566400    0.935277939    2.108163834
        1 ACE   HH32       3    1.074337721    1.037187934    2.206760406
        1 ACE   HH33       4    1.048208833    1.004362464    2.024867773
        1 ACE   C          5    1.194072366    1.147308826    2.073846102
        1 ACE   O          6    1.142704129    1.229152918    1.999127984
        2 ALA   N          7    1.309494257    1.162267447    2.132897854
        2 ALA   HN         8    1.339519024    1.092915297    2.196325302
        2 ALA   CA         9    1.403393507    1.269239664    2.110365152
        2 ALA   HA        10    1.377287865    1.325353146    2.023513317
        2 ALA   CB        11    1.414769530    1.358132243    2.234688044
        2 ALA   HB1       12    1.440894365    1.293333054    2.326516151
        2 ALA   HB2       13    1.493342519    1.434360743    2.223044872
        2 ALA   HB3       14    1.315252423    1.413229465    2.254475832
        2 ALA   C         15    1.537224412    1.212282062    2.064972878
        2 ALA   O         16    1.580033302    1.103038311    2.096084595
        3 ARG   N         17    1.617468953    1.286837935    1.985834837
        3 ARG   HN        18    1.577086687    1.382247448    1.962693930
        3 ARG   CA        19    1.753297091    1.257931352    1.948590279
        ...

Se visualiza `arv.log`.

    Initial temperature: 298.011 K

    Started mdrun on rank 0 Mon Mar 10 17:33:20 2025
               Step           Time
                  0        0.00000

       Energies (kJ/mol)
               Bond          Angle            U-B    Proper Dih.  Improper Dih.
        3.69924e+03    2.02073e+03    1.14815e+02    5.68654e+01    1.04274e+01
          CMAP Dih.          LJ-14     Coulomb-14        LJ (SR)  Disper. corr.
       -8.58925e+00    2.92018e+01   -1.05816e+03    7.15782e+03   -1.55883e+02
       Coulomb (SR)   Coul. recip.      Potential    Kinetic En.   Total Energy
       -4.81137e+04    1.44763e+02   -3.61024e+04    9.82123e+03   -2.62812e+04
      Conserved En.    Temperature Pres. DC (bar) Pressure (bar)
       -2.62812e+04    2.99082e+02   -9.58703e+01    1.91375e+02

               Step           Time
                 20        0.01000

      Energies (kJ/mol)
               Bond          Angle            U-B    Proper Dih.  Improper Dih.
        3.68335e+03    1.96341e+03    1.26291e+02    6.06098e+01    1.22199e+01
          CMAP Dih.          LJ-14     Coulomb-14        LJ (SR)  Disper. corr.
       -4.54557e+00    2.30252e+01   -1.06894e+03    7.06760e+03   -1.55883e+02
       Coulomb (SR)   Coul. recip.      Potential    Kinetic En.   Total Energy
       -4.79755e+04    1.44451e+02   -3.61239e+04    9.90597e+03   -2.62180e+04
      Conserved En.    Temperature Pres. DC (bar) Pressure (bar)
       -2.62817e+04    3.01662e+02   -9.58703e+01    1.45345e+02

               Step           Time
                 40        0.02000
    ...

La simulación de **equilibración NVT** en **GROMACS** se configura con
**1,000,000 pasos** (500 ps), asegurando el mantenimiento de la
temperatura en **298K**. Se copian los archivos esenciales (`arv.top`,
`arv.g96`, `runNVT.mdp`), estableciendo la frecuencia de guardado de
coordenadas, velocidades y energías cada **10 fs**. La salida de
ejecución (`slurm-43314.out`) confirma que la simulación finalizó
correctamente con un rendimiento de **22.561 ns/día**, y la inspección
del archivo `arv.g96` muestra la disposición final de los átomos en el
sistema. El análisis de `arv.log` refleja valores de **energía potencial
estable en torno a -36,100 kJ/mol**, con una **energía cinética de
~9,900 kJ/mol** y una temperatura oscilando alrededor del valor
esperado. Además, se observa una presión fluctuante sin tendencias
extremas, lo que indica que el sistema alcanzó un **estado de equilibrio
térmico**, listo para iniciar la simulación de producción.

### **3.6. Tarea extra 2:**

**Extraer los datos de cada residuo a 298K y 400K**

    grep "ALA-2" rama_arv.xvg | awk '{print $1, $2}' > rama_ALA2.dat
    grep "ARG-3" rama_arv.xvg | awk '{print $1, $2}' > rama_ARG3.dat
    grep "VAL-4" rama_arv.xvg | awk '{print $1, $2}' > rama_VAL4.dat

    awk '{printf "%.1f %.1f\n", $1, $2}' rama_ALA2.dat | sort -n | uniq -c | awk '{print $2, $3, $1}' > rama_density_ALA2.dat
    awk '{printf "%.1f %.1f\n", $1, $2}' rama_ARG3.dat | sort -n | uniq -c | awk '{print $2, $3, $1}' > rama_density_ARG3.dat
    awk '{printf "%.1f %.1f\n", $1, $2}' rama_VAL4.dat | sort -n | uniq -c | awk '{print $2, $3, $1}' > rama_density_VAL4.dat

    awk '{x=int(($1+180)/10); y=int(($2+180)/10); grid[x,y]++} 
         END {for (i in grid) {split(i, coords, SUBSEP); print coords[1]*10-180, coords[2]*10-180, grid[i]}}' rama_arv.xvg > rama_density.dat

    reset

    # Configuración general
    set terminal pngcairo size 800,600

    # Configuración de la paleta de colores (amarillo → blanco → violeta)
    set palette defined (0 "yellow", 0.5 "white", 1 "violet")

    # Ajustar rango de color para que vaya de negativo a 0
    set cbrange [-0.002:0]

    # Configuración de los ejes
    set xlabel "Φ (Phi) [°]"
    set ylabel "Ψ (Psi) [°]"
    set xrange [-180:180]
    set yrange [-180:180]
    set grid

    # Ajuste de interpolación para mejorar la visualización
    set pm3d interpolate 2,2

    # Tamaño de los puntos para mejorar la visualización de la nube
    set pointsize 1.5
    set style data points

    ############# PLOT PARA ALA-2 #############
    set output 'ramachandran_ALA2_nube.png'
    set title "Diagrama de Ramachandran - ALA-2"
    plot "rama_density_ALA2.dat" using 1:2:3 with points pt 7 lc palette notitle

    ############# PLOT PARA ARG-3 #############
    set output 'ramachandran_ARG3_nube.png'
    set title "Diagrama de Ramachandran - ARG-3"
    plot "rama_density_ARG3.dat" using 1:2:3 with points pt 7 lc palette notitle

    ############# PLOT PARA VAL-4 #############
    set output 'ramachandran_VAL4_nube.png'
    set title "Diagrama de Ramachandran - VAL-4"
    plot "rama_density_VAL4.dat" using 1:2:3 with points pt 7 lc palette notitle

El procedimiento realizado tiene como objetivo extraer y visualizar los
datos de los ángulos Phi (Φ) y Psi (Ψ) de cada residuo en la
**simulación a 298K y 400K**, generando los diagramas de Ramachandran
para **los residuos ALA-2, ARG-3 y VAL-4**. Primero, se extraen los
valores de Φ y Ψ de cada residuo del archivo `rama_arv.xvg` mediante
`grep`, guardándolos en archivos individuales (`rama_ALA2.dat`,
`rama_ARG3.dat`, `rama_VAL4.dat`). Luego, se procesan estos datos para
calcular la densidad de puntos en el espacio Ramachandran mediante
`awk`, agrupando coordenadas y contando su frecuencia
(`rama_density_ALA2.dat`, `rama_density_ARG3.dat`,
`rama_density_VAL4.dat`). Adicionalmente, se genera una malla de
densidad global (`rama_density.dat`) para representar la distribución de
conformaciones. Finalmente, se configuran y generan los gráficos con
`gnuplot`, estableciendo un esquema de colores que va de amarillo a
violeta, interpolación de datos para mejorar la visualización y
representación de la densidad de puntos como una nube. Se generan tres
diagramas de Ramachandran (`ramachandran_ALA2_nube.png`,
`ramachandran_ARG3_nube.png`, `ramachandran_VAL4_nube.png`), permitiendo
comparar la distribución conformacional de los residuos en ambas
temperaturas.

**- A 298K.**

<figure>
<img
src="/home/alumno08/MM/tarea/4-analysis/280k/ramachandran_ALA2_nube.png"
alt="Diagrama de Ramachandran - ALA-2 298K" />
<figcaption aria-hidden="true">Diagrama de Ramachandran - ALA-2
298K</figcaption>
</figure>

<figure>
<img
src="/home/alumno08/MM/tarea/4-analysis/280k/ramachandran_ARG3_nube.png"
alt="Diagrama de Ramachandran - ARG-3 298K" />
<figcaption aria-hidden="true">Diagrama de Ramachandran - ARG-3
298K</figcaption>
</figure>

<figure>
<img
src="/home/alumno08/MM/tarea/4-analysis/280k/ramachandran_VAL4_nube.png"
alt="Diagrama de Ramachandran - VAL-4 298K" />
<figcaption aria-hidden="true">Diagrama de Ramachandran - VAL-4
298K</figcaption>
</figure>

**- A 400K**

<figure>
<img
src="/home/alumno08/MM/tarea/4-analysis/400k/ramachandran_ALA2_nube.png"
alt="Diagrama de Ramachandran - ALA-2 400K" />
<figcaption aria-hidden="true">Diagrama de Ramachandran - ALA-2
400K</figcaption>
</figure>

![Diagrama de Ramachandran - ARG-3
400K](/home/alumno08/MM/tarea/4-analysis/400k/ramachandran_ARG3_nube.png)
![Diagrama de Ramachandran - VAL-4
400K](/home/alumno08/MM/tarea/4-analysis/400k/ramachandran_VAL4_nube.png)

**Análisis de los resultados.**

La comparación de los **diagramas de Ramachandran** entre **298K**
(primeras tres imágenes) y **400K** (últimas tres imágenes) revela
diferencias en la distribución conformacional de los residuos **ALA-2,
ARG-3 y VAL-4**, reflejando el impacto del aumento de temperatura en la
flexibilidad estructural. A **298K**, las conformaciones de **ALA-2 y
ARG-3** se encuentran más concentradas en regiones definidas, mientras
que **VAL-4** muestra una distribución compacta en torno a valores
negativos de **Φ y Ψ**. A **400K**, se observa una mayor dispersión en
todas las estructuras, lo que sugiere un aumento en la **flexibilidad
conformacional** debido a la mayor energía térmica. En particular, los
residuos tienden a explorar un espacio conformacional más amplio, lo que
podría estar asociado con una menor estabilidad estructural y una mayor
propensión a fluctuaciones dinámicas.

También se realiza un **análisis estadístico de los ángulos Phi (ϕ) y
Psi (ψ)** de residuos de una proteína en diferentes temperaturas (298K y
400K) a partir de diagramas de Ramachandran.

    # Definir rutas de los archivos
    files_298K <- list(
      "ALA-2" = "/home/alumno08/MM/tarea/4-analysis/280k/rama_ALA2.dat",
      "ARG-3" = "/home/alumno08/MM/tarea/4-analysis/280k/rama_ARG3.dat",
      "VAL-4" = "/home/alumno08/MM/tarea/4-analysis/280k/rama_VAL4.dat"
    )

    files_400K <- list(
      "ALA-2" = "/home/alumno08/MM/tarea/4-analysis/400k/rama_ALA2.dat",
      "ARG-3" = "/home/alumno08/MM/tarea/4-analysis/400k/rama_ARG3.dat",
      "VAL-4" = "/home/alumno08/MM/tarea/4-analysis/400k/rama_VAL4.dat"
    )

    # Función para calcular estadísticas angulares
    calcular_estadisticas <- function(df) {
      data.frame(
        "Phi Range" = max(df$Phi) - min(df$Phi),
        "Psi Range" = max(df$Psi) - min(df$Psi),
        "Phi Mean" = mean(df$Phi),
        "Psi Mean" = mean(df$Psi),
        "Phi SD" = sd(df$Phi),
        "Psi SD" = sd(df$Psi)
      )
    }

    # Leer los archivos y calcular estadísticas para 298K
    estadisticas_298K <- lapply(files_298K, function(filepath) {
      df <- read.table(filepath, header=FALSE, col.names=c("Phi", "Psi"))
      calcular_estadisticas(df)
    })

    # Leer los archivos y calcular estadísticas para 400K
    estadisticas_400K <- lapply(files_400K, function(filepath) {
      df <- read.table(filepath, header=FALSE, col.names=c("Phi", "Psi"))
      calcular_estadisticas(df)
    })

    # Convertir a dataframes
    df_estadisticas_298K <- do.call(rbind, estadisticas_298K)
    df_estadisticas_400K <- do.call(rbind, estadisticas_400K)

    # Agregar nombres de fila
    rownames(df_estadisticas_298K) <- names(files_298K)
    rownames(df_estadisticas_400K) <- names(files_400K)

    # Mostrar resultados
    print("Estadísticas de Ángulos a 298K")

    ## [1] "Estadísticas de Ángulos a 298K"

    print(df_estadisticas_298K)

    ##       Phi.Range Psi.Range   Phi.Mean  Psi.Mean    Phi.SD    Psi.SD
    ## ALA-2   71.6542   36.6640  -66.88251 152.32714 17.279790  7.213724
    ## ARG-3   49.0044   35.9790  -65.71950 135.81348  7.559887  7.288491
    ## VAL-4   40.5738   47.0776 -113.91755 -36.53358  9.285001 10.388131

    print("Estadísticas de Ángulos a 400K")

    ## [1] "Estadísticas de Ángulos a 400K"

    print(df_estadisticas_400K)

    ##       Phi.Range Psi.Range   Phi.Mean Psi.Mean    Phi.SD    Psi.SD
    ## ALA-2   70.2887   60.5800   73.74990 152.8121 12.932755 14.988433
    ## ARG-3   48.5860   46.7130 -148.97357 160.6456 10.267817  8.113099
    ## VAL-4   47.9193   60.0976  -87.87824 120.5674  9.679041 13.328673

A **400K**, el tripéptido ARV muestra una mayor **flexibilidad
conformacional** en comparación con **298K**, reflejada en el aumento
del **rango y la desviación estándar** de los ángulos **Psi** en todos
los residuos, especialmente en **ALA-2 y VAL-4**. Además, los valores
medios de **Phi y Psi** presentan cambios drásticos, con **ARG-3**
mostrando una alteración significativa en **Phi.Mean** (de -65.72 a
-148.97) y **VAL-4** en **Psi.Mean** (de -36.53 a 120.57), indicando una
reorientación estructural importante. La mayor dispersión angular y los
cambios en la media sugieren que el tripéptido **explora un espacio
conformacional más amplio a temperaturas elevadas**, lo que podría
afectar su estabilidad y función biológica en condiciones extremas.

## **4. Conclusiones finales**

El presente estudio ha permitido caracterizar el comportamiento dinámico
y conformacional del **tripéptido ARV (Alanina-Arginina-Valina)** a
través de simulaciones de **dinámica molecular** a diferentes
temperaturas (298K y 400K). A partir de los análisis realizados, se han
obtenido los siguientes puntos clave:

### **4.1. Estabilidad estructural y efecto de la temperatura:**

-   A **298K**, el tripéptido se mantiene relativamente estable, con
    fluctuaciones limitadas en la distancia entre átomos clave (CA-CB) y
    en los radios de giro.

-   A **400K**, se observan **fluctuaciones más pronunciadas**, lo que
    sugiere un aumento en la flexibilidad estructural y posibles
    indicios de desestabilización térmica.

-   El análisis de los **diagramas de Ramachandran** evidencia una mayor
    dispersión conformacional a 400K, indicando una mayor exploración
    del espacio conformacional del péptido debido al incremento en la
    energía térmica.

### **4.2. Propiedades dinámicas y flexibilidad conformacional:**

-   Se ha observado que **la temperatura afecta la estabilidad de los
    ángulos Phi y Psi**, con mayores variaciones a 400K, lo que sugiere
    una mayor libertad conformacional.

-   Los **ángulos dihédricos** presentan oscilaciones más amplias a
    temperaturas elevadas, reforzando la idea de que el tripéptido
    adquiere mayor flexibilidad estructural con el aumento de
    temperatura.

### **4.3. Interacciones moleculares y puentes de hidrógeno:**

-   A 298K, se detectaron **puentes de hidrógeno intramoleculares**
    estables, contribuyendo a la conservación de la estructura del
    tripéptido.

-   A 400K, estos **puentes de hidrógeno se reducen notablemente**,
    indicando una posible pérdida de estabilidad estructural a
    temperaturas más altas.

### **4.4. Energía y dinámica del sistema:**

-   La **energía total** del sistema se mantiene relativamente estable
    en ambas condiciones, aunque **a 400K se observa una mayor
    variabilidad**, lo que sugiere un incremento en la movilidad
    molecular.

-   La **energía cinética** es mayor a temperaturas elevadas, lo que es
    consistente con el **mayor grado de agitación térmica** observado en
    las simulaciones.
