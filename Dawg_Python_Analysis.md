# A Brief Analysis of Dawg as a Python Module
 - The goal of these changes is to use DAWG effectively in a high-level language like Python
 - Use DAWG in Python should make it easier to write tests and integrate with other Python libraries
 - Boost.Python, Cython, and SWIG are the three main wrapping tools analyzed
 - A BASH script was used to do the time trials **analyze_wrappers.sh**
    - Each wrapper was tested with the **basic-dna.dawg** trick file
    - Each wrapper used 10 reps and seed **212121**
    - Each wrapper program was called continuously 10000 times, output redirected to **/dev/null**
    - Linux **time** command was used to record the time (real time)
 - The sizes of their libraries were analyzed with **ls -lh**
 - A simple API was created to try to show how to use DAWG quickly from Python
    - This API was more or less consistent between all three wrappers
    - It consisted of a **run** command which took a DAWG trick file, output, rep, and seed
    - There are some discrepencies between the APIs and how they are called:
        - For instance, Cython requires that a string be converted into bytes before getting called
 - For the API example, we used the Cython generated wrapper

## Boost.Python
 - Time trial: 6min45.847s
 - **dawg.so** size: 6.7M

## Cython
 - Time trial: 3m51.320s
 - **PyDawg.cpython-36m-x86_64-linux-gnu.so** size: 20M

## SWIG
 - Time trial: 7m17.167
 - **_dawg.so** size: 6.6M
