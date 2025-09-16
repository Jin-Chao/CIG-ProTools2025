#/bin/bash

if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 <CIG_DIR> <GENERATOR_DIR>"
    exit 1
fi

CIG_DIR=$1
GENERATOR_DIR="$2"

"${GENERATOR_DIR}"/cig_generator.in --cig ${CIG_DIR}/cacheusage.cig.conflict --function initialization --output Figure_4.pdf
"${GENERATOR_DIR}"/cig_generator.in --cig ${CIG_DIR}/cacheusage.cig.conflict --function conflict --output Figure_5.pdf

"${GENERATOR_DIR}"/cig_generator.in --cig ${CIG_DIR}/cacheusage.cig.matrix --function simple --output Figure_6_L.pdf
"${GENERATOR_DIR}"/cig_generator.in --cig ${CIG_DIR}/cacheusage.cig.matrix --function simple --threshold 0.01 --output Figure_6_R.pdf
"${GENERATOR_DIR}"/cig_generator.in --cig ${CIG_DIR}/cacheusage.cig.matrix --function blocked --threshold 0.03 --output Figure_7.pdf

"${GENERATOR_DIR}"/cig_generator.in --cig ${CIG_DIR}/cacheusage.cig.bmt --function jacobi --miss Conflict --output Figure_8.pdf
"${GENERATOR_DIR}"/cig_generator.in --cig ${CIG_DIR}/cacheusage.cig.bmt_pad --function jacobi --threshold 0.01 --output Figure_9.pdf
"${GENERATOR_DIR}"/cig_generator.in --cig ${CIG_DIR}/cacheusage.cig.bmt_aos --function jacobi --threshold 0.01 --output Figure_10.pdf

"${GENERATOR_DIR}"/cig_generator.in --cig ${CIG_DIR}/cacheusage.cig.irsmk --function rmat --miss Conflict --threshold 0.01 --output Figure_12.pdf
"${GENERATOR_DIR}"/cig_generator.in --cig ${CIG_DIR}/cacheusage.cig.irsmk_aos --function rmat  --threshold 0.01 --output Figure_13.pdf
"${GENERATOR_DIR}"/cig_generator.in --cig ${CIG_DIR}/cacheusage.cig.correlation --function correlation  --threshold 0.01 --output Figure_14.pdf
