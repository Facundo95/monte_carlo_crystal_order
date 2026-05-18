#!/bin/bash

# Comprobamos que se pasen al menos los 2 argumentos necesarios
if [ "$#" -lt 2 ]; then
    echo "Uso: $0 <elemento_aleatorio> <archivo.xyz>"
    echo "Ejemplo: $0 Mn sistema.xyz"
    exit 1
fi

target_element="$1"
input_file="$2"
output_file="procesado_${input_file}"

awk -v target="$target_element" '
BEGIN { srand() } 

{
    # Si la fila es menor a 3, imprimir tal cual
    if (NR < 3) {
        print $0
    } 
    else {
        # Por defecto el valor es 0.0
        val = "0.0"
        
        # Si el primer campo coincide con el elemento pasado por comando
        if ($1 == target) {
            val = (rand() < 0.5) ? "1.0" : "-1.0"
        }
        
        # Imprimir la línea original seguida de la nueva columna
        print $0, val
    }
}' "$input_file" > "$output_file"

echo "Procesado elemento: $target"
echo "Archivo guardado como: $output_file"
