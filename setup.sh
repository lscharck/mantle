#!/bin/bash

mkdir -p lib obj

sed -i '' -e 's/#define REB_RESTRICT restrict/#define REB_RESTRICT __restrict__/g' ./rebound/src/rebound.h

