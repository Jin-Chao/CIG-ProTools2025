EXP_DIR=../..
PROJ_DIR=${EXP_DIR}/Simple

SRC_DIR=${PROJ_DIR}/src
BIN_DIR=${PROJ_DIR}/bin
OBJ_DIR=${PROJ_DIR}/obj
UTIL_DIR=${EXP_DIR}/Util

INC=-I${UTIL_DIR}

CC=gcc
CFLAGS=-O3 -g -Wall
DEPS=${UTIL_DIR}/util.h ${UTIL_DIR}/util.c

${OBJ_DIR}/%.o: %.c 
	${CC} ${CFLAGS} ${INC} -c -o $@ $< 

${OBJ_DIR}/%.o: ${UTIL_DIR}/%.c 
	${CC} ${CFLAGS} ${INC} -c -o $@ $<

all: matrix conflict

matrix: ${OBJ_DIR}/matrix.o ${OBJ_DIR}/util.o
	${CC} ${CFLAGS} -o ${BIN_DIR}/$@ $^ 

conflict: ${OBJ_DIR}/conflict.o ${OBJ_DIR}/util.o
	${CC} ${CFLAGS} -o ${BIN_DIR}/$@ $^ 

.PHONY: clean

clean:
	rm ${BIN_DIR}/* ${OBJ_DIR}/*

echos:
	@echo "INC: " ${INC}
