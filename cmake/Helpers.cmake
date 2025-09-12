include_guard()

function(set_output_path lib path)
    set_target_properties(${lib} PROPERTIES
	PREFIX ""
	SUFFIX "${PYTHON_MODULE_SUFFIX}"
	OUTPUT_NAME "${lib}"
        LIBRARY_OUTPUT_DIRECTORY "${path}"
        RUNTIME_OUTPUT_DIRECTORY "${path}"
    )

    foreach(CONFIG IN ITEMS RELEASE DEBUG RELWITHDEBINFO MINSIZEREL)
        set_target_properties(${lib} PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY_${CONFIG} "${path}"
            RUNTIME_OUTPUT_DIRECTORY_${CONFIG} "${path}"
        )
    endforeach()
endfunction()

