#!/usr/bin/env bash

set -euo pipefail

usage() {
	cat <<'EOF'
Uso: from_fl_to_xyz.sh <archivo.fl> [archivo.xyz] [Atom1] [Atom2] [Atom3]

Convierte un archivo .fl a un archivo .xyz usando la misma decodificacion que
el codigo C++ del proyecto.

Cada valor .fl se transforma en un par (especie, spin) y se escribe como:
	atom x y z spin

Las especies se mapean con Atom1, Atom2 y Atom3.

Por defecto el archivo de salida se genera cambiando la extension a .xyz.
Si no se indican nombres atomicos, se usan Atom1, Atom2 y Atom3.
EOF
}

if [[ $# -lt 1 || $# -gt 5 ]]; then
	usage
	exit 1
fi

input_file="$1"
output_file="${2:-${input_file%.*}.xyz}"
atom1="${3:-Atom1}"
atom2="${4:-Atom2}"
atom3="${5:-Atom3}"

if [[ ! -f "$input_file" ]]; then
	echo "Error: no existe el archivo de entrada: $input_file" >&2
	exit 1
fi

awk -v atom1="$atom1" -v atom2="$atom2" -v atom3="$atom3" '
function trim(s) {
	gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
	return s
}

{
	line = trim($0)
	if (line == "" || line ~ /^#/) {
		next
	}

	if (line !~ /^-?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$/) {
		printf "Error: linea no numerica en .fl: %s\n", $0 > "/dev/stderr"
		exit 1
	}

	values[++count] = line + 0
}

END {
	if (count == 0) {
		print "Error: el archivo .fl no contiene datos" > "/dev/stderr"
		exit 1
	}

	side = int(exp(log(count / 2.0) / 3.0) + 0.5)
	if (2 * side * side * side != count) {
		printf "Error: cantidad de sitios incompatible con una red side x side x (2*side): %d\n", count > "/dev/stderr"
		exit 1
	}

	print count
	print ""

	plane_size = side * (2 * side)

	for (i = 1; i <= count; ++i) {
		value = values[i]

		spin = int((-1.0 / 12.0) * value * value * value - (1.0 / 3.0) * value * value + (7.0 / 12.0) * value + (5.0 / 6.0))
		species = int((-1.0 / 12.0) * value * value * value + (1.0 / 3.0) * value * value + (7.0 / 12.0) * value - (5.0 / 6.0))

		if (species == 1) {
			atom = atom1
		} else if (species == 0) {
			atom = atom2
		} else if (species == -1) {
			atom = atom3
		} else {
			printf "Error: estado .fl no reconocido en la posicion %d: valor=%g especie=%d spin=%d\n", i, value, species, spin > "/dev/stderr"
			exit 1
		}

		site = i - 1
		x = int(site / plane_size)
		tmp = site - (x * plane_size)
		y = int(tmp / (2 * side))
		z = tmp - (y * (2 * side))

		if ((z % 2) != 0) {
			printf "%s %.1f %.1f %.1f %d\n", atom, x + 0.5, y + 0.5, z / 2.0, spin
		} else {
			printf "%s %.1f %.1f %.1f %d\n", atom, x, y, z / 2.0, spin
		}
	}
}
' "$input_file" > "$output_file"

echo "Archivo convertido: $output_file"
