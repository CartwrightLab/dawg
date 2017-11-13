PYTHON_ADD_MODULE(dawg_python dawg.cpp)
FILE(COPY dawg.py DESTINATION .)
ADD_TEST(NAME boost.python COMMAND ${PYTHON_EXECUTABLE} dawg.py)
