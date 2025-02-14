#!/bin/bash

mkdir -pv lib obj

echo "applying patch"
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
	sed -i 's/#define REB_RESTRICT restrict/#define REB_RESTRICT __restrict__/g' ./rebound/src/rebound.h

	sed -i 's/LIB+= -lm -lrt/LIB+= -lm -lrt -lpthread/g' ./rebound/src/Makefile.defs
elif [[ "$OSTYPE" == "darwin"* ]]; then
	sed -i '' -e 's/#define REB_RESTRICT restrict/#define REB_RESTRICT __restrict__/g' ./rebound/src/rebound.h

	sed -i '' -e 's/LIB+= -lm -lrt/LIB+= -lm -lrt -lpthread/g' ./rebound/src/Makefile.defs
fi

