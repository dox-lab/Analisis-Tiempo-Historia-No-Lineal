def ReadSMDFile(inFilename, outFilename, dt):
    # Abre el archivo de entrada y captura el error si no puede ser leído
    try:
        inFileID = open(inFilename, 'r')
    except IOError:
        print(f"Cannot open {inFilename} for reading", file=sys.stderr)
        return

    # Abre el archivo de salida para escritura
    outFileID = open(outFilename, 'w')

    # Indicador que señala cuando se ha encontrado el dt y se deben leer
    # los valores de movimiento de suelo -- ASUME que dt está en la última
    # línea del encabezado
    flag = 0

    # Lee cada línea en el archivo
    for line in inFileID:
        line = line.strip()  # Elimina espacios en blanco al inicio y final de la línea

        if len(line) == 0:
            # Línea en blanco --> no hacer nada
            continue
        elif flag == 1:
            # Imprime los valores de movimiento de suelo en el archivo de salida
            outFileID.write(line + '\n')
        else:
            # Buscar el dt en las líneas del encabezado
            for word in line.split():
                # Leer el paso de tiempo
                if flag == 1:
                    dt = word
                    break
                # Encontrar el token deseado y establecer el indicador
                if word.startswith("DT="):
                    flag = 1
    # Cierra el archivo de salida
    outFileID.close()

    # Cierra el archivo de entrada
    inFileID.close()