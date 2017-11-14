# Ignoring FIND_PACKAGE checks for now
EXECUTE_PROCESS(COMMAND "sh swigify.sh")
FILE(COPY example.py _example.so DESTINATION .)
ADD_TEST(NAME swig COMMAND ${PYTHON_EXECUTABLE} example.py)
