EXECUTE_PROCESS(COMMAND "${PYTHON_EXECUTABLE} setup.py build_ext -if")
FILE(COPY DawgDriver.py PyDawg.cpp PyDawg.cpython-36m-x86_64-linux-gnu.so DESTINATION .)
ADD_TEST(NAME cython COMMAND ${PYTHON_EXECUTABLE} DawgDriver.py)
