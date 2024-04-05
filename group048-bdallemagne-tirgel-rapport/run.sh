#!/bin/bash
ORDER="ProjectPreProcessor Project ProjectPostProcessor"
GREEN='\033[0;32m'
NC='\033[0m' # No Color
RED='\033[0;31m'
#usr/bin/cmake --no-warn-unused-cli -DCMAKE_BUILD_TYPE:STRING=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=TRUE -DCMAKE_C_COMPILER:FILEPATH=/usr/bin/gcc -DCMAKE_CXX_COMPILER:FILEPATH=/usr/bin/g++ -S/elem_fini/LEPL1110_elem_fini/group048-bdallemagne-tirgel-rapport/Project/Project -B/elem_fini/LEPL1110_elem_fini/group048-bdallemagne-tirgel-rapport/Project/Project/build -G "Unix Makefiles"
MAINDIR=$(pwd)

echo -e "$GREEN [LOG]:$NC this will compile and run the 3 main components of the project"

# Iterate over ORDER
for i in $ORDER ; do
    if [ $i == "Project" ] || [ $i == "ProjectPostProcessor" ]; then 
        # Copying the data in the next folder
        echo -e "$GREEN [LOG]:$NC copying data to $i"
        cp -r ../data ../../$i
        cd ../..
    fi
    echo -e "$GREEN [LOG]:$NC compiling $i"
    cd $i
    cmake --no-warn-unused-cli -DCMAKE_BUILD_TYPE:STRING=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=TRUE -DCMAKE_C_COMPILER:FILEPATH=/usr/bin/gcc -DCMAKE_CXX_COMPILER:FILEPATH=/usr/bin/g++ -S$MAINDIR/$i -B$MAINDIR/$i/build -G "Unix Makefiles" > logs.txt
    # Check succes of last command
    if [ $? -ne 0 ]; then
        echo -e "$RED [ERROR]:$NC cmake failed for $i, please check the logs"
        exit 1
    fi
    cmake --build $MAINDIR/$i/build --config Debug --target all -j 18 -- > logs.txt
    cd build
    ./myFem
done

echo -e "$GREEN [LOG]:$NC Finished"