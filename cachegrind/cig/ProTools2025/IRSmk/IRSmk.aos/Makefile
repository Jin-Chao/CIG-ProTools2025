EXP_DIR=../..
PROJ_DIR=${EXP_DIR}/IRSmk
  
SRC_DIR=${PROJ_DIR}/IRSmk.aos
BIN_DIR=${PROJ_DIR}/bin
OBJ_DIR=${PROJ_DIR}/obj
UTIL_DIR=${EXP_DIR}/Util

INC=-I${UTIL_DIR}

CC = gcc 
IRSmk_PROG = ${BIN_DIR}/irsmk.aos
IRSmk_OBJS = ${OBJ_DIR}/main.o ${OBJ_DIR}/rmatmult3.o ${OBJ_DIR}/utility.o ${OBJ_DIR}/util.o
CFLAGS = -fopenmp -O3 -g
LDFLAGS = -lm

RM = rm

${OBJ_DIR}/%.o: %.c 
	${CC} ${CFLAGS} ${INC} -c -o $@ $<

${OBJ_DIR}/%.o: ${UTIL_DIR}/%.c 
	${CC} ${CFLAGS} ${INC} -c -o $@ $<

all : ${IRSmk_PROG}

${IRSmk_PROG} : ${IRSmk_OBJS}
	${CC} ${CFLAGS} -o $@ ${IRSmk_OBJS} ${LDFAGS}
	${RM} -f ${IRSmk_OBJS} 

clean :
	${RM} -f ${IRSmk_PROG} ${IRSmk_OBJS}
