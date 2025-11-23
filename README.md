# Análisis Tiempo-Historia No Lineal de un Edificio de 8 Pisos

## Descripción

Este proyecto realiza el análisis tiempo-historia no lineal de una estructura de concreto armado de 8 pisos usando programación interactiva con **Python** y la librería **OpenSeesPy**. El análisis se basa en las directrices de la norma técnica **E.030** de Perú para sismos y compara los resultados con el software comercial **ETABS**. Además, se utiliza el formato **DXF** para importar el diseño estructural inicial y **EZDXF** para manipular estos archivos en Python.

## Objetivos

- Implementar un análisis estructural no lineal en tiempo usando Python y OpenSeesPy.
- Modelar una estructura de concreto armado irregular con muros usando archivos **DXF**.
- Realizar un análisis sísmico siguiendo las directrices de la norma **E.030**.
- Comparar los resultados del análisis con **ETABS**.
- Explorar las capacidades de OpenSeesPy para modelado y simulación sísmica avanzada.

## Estructura del Proyecto

El proyecto está organizado en las siguientes carpetas:

```plaintext
📁 tesis-analisis-tiempo-historia-no-lineal
│
├── 📁 data                     # Archivos de entrada como DXF, acelerogramas y propiedades de materiales
├── 📁 etabs                    # Modelos y resultados obtenidos con ETABS
├── 📁 src                      # Código fuente organizado en modelos, análisis y comparación
├── 📁 docs                     # Documentos como la tesis en PDF y presentaciones
├── 📁 resultados               # Gráficas y resultados exportados
└── README.md                   # Este archivo
```

## Requisitos

Las siguientes herramientas y bibliotecas son necesarias para ejecutar el proyecto:

- **Python 3.11**
- **OpenSeesPy**
- **EZDXF**
- **NumPy**
- **Matplotlib**

Puedes instalar las dependencias ejecutando el siguiente comando:

```bash
pip install -r requirements.txt
```

## Uso del Proyecto

### Clonar el repositorio

Primero, clona el repositorio en tu máquina local:

```bash
git clone https://github.com/TU_USUARIO/analisis-tiempo-historia-no-lineal.git
cd analisis-tiempo-historia-no-lineal
```

### Activar el entorno virtual

Activa el entorno de conda (o crea uno si no lo tienes):

```bash
conda activate analisis_no_lineal_env
```

### Ejecutar el análisis

Puedes ejecutar el análisis no lineal en tiempo corriendo el script principal en el directorio `src`:

```bash
python src/analisis/analisis_sismico.py
```

### Comparar resultados con ETABS

Para comparar los resultados obtenidos con OpenSeesPy y ETABS:

```bash
python src/comparacion/comparacion_etabs.py
```

## Documentación

La documentación completa del proyecto, incluyendo la metodología y los resultados, se encuentra en el archivo `docs/tesis.pdf`.

## Contacto

Para preguntas o sugerencias, puedes contactarme en [daniel.medina@utec.edu.pe].

## Licencia

Este proyecto está bajo la licencia **MIT**. Consulta el archivo [LICENSE](LICENSE) para más detalles.
