# Invoked by add_custom_command to embed a binary/text file as a C byte array.
#
# Required variables (passed via -D):
#   INPUT   - path to the source file
#   OUTPUT  - path to write the generated C include
#   VARNAME - C variable name to use

if(NOT EXISTS "${INPUT}")
    message(FATAL_ERROR "Cannot embed missing file: ${INPUT}")
endif()

# read the file as one long string of hexadecimal digits, two per byte
file(READ "${INPUT}" hex HEX)

string(LENGTH "${hex}" hex_length)
math(EXPR num_bytes "${hex_length} / 2")

# turn every pair of digits into a C byte literal: "2f2f20" -> "0x2f,0x2f,0x20,"
string(REGEX REPLACE "([0-9a-f][0-9a-f])" "0x\\1," body "${hex}")

# drop the comma trailing the last byte
string(REGEX REPLACE ",$" "" body "${body}")

file(WRITE "${OUTPUT}"
    "unsigned char ${VARNAME}[] = { ${body} };\n"
    "unsigned int ${VARNAME}_len = ${num_bytes};\n"
)
