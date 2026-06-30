# Invoked by add_custom_command to embed a binary/text file as a C byte array.
# Required variables (passed via -D):
#   INPUT   - path to the source file
#   OUTPUT  - path to write the generated C include
#   VARNAME - C variable name to use (overrides xxd's path-derived name)

execute_process(
    COMMAND xxd -i -n "${VARNAME}" "${INPUT}"
    OUTPUT_FILE "${OUTPUT}"
    RESULT_VARIABLE result
)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "xxd failed (exit code ${result}) for input: ${INPUT}")
endif()
