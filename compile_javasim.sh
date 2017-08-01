#!/bin/sh

mkdir -p reference_generation/bin
lib=reference_generation/lib/
javac -cp ${lib}commons-lang3-3.1.jar:${lib}/picard-1.77.jar:lib/sam-1.77.jar -d reference_generation/bin reference_generation/src/*.java
