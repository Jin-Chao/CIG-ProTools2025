# README: CIG Example Suite created for SC ProTools 2025

This README describes how to build and run the examples in the Cache Interaction Graph (CIG) suite. A master Makefile in the `examples/` directory builds all sub-projects.

---

## 1. Simple Examples

**Directory:** `Simple/`

* `conflict`: Demonstrates conflict cache misses by mapping three variables to the same cache set.
* `matrix`: Compares naive and blocked matrix multiplication.

**Usage:**

```bash
cd Simple/src
make

# Run examples
../bin/conflict -s 1024
../bin/matrix -s 256 -b 8
../bin/matrix -s 512 -b 8
../bin/matrix -s 1024 -b 8
```

---

## 2. Himeno Benchmark

**Source:** [RIKEN HimenoBMT](https://i.riken.jp/en/supercom/documents/himenobmt/download/win-mac/)

* `himenoBMTxpa.c`: Original version with memory tracking
* `himenoBMTxpa_pad.c`: Version with address shifting
* `himenoBMTxpa_aos.c`: Optimized using Array of Structures (AoS)

**Usage:**

```bash
cd Himeno/src
make   # Builds all versions with both O0 and O3 optimizations

# Run examples
../bin/bmt.O0 -s M -l 100
../bin/bmt_pad.O0 -s M -l 100
../bin/bmt_aos.O0 -s M -l 100
../bin/bmt.O3 -s M -l 100
../bin/bmt_pad.O3 -s M -l 100
../bin/bmt_aos.O3 -s M -l 100
```

---

## 3. IRSmk Benchmark

**Source:** [Microsoft IRSmk Benchmark](https://github.com/microsoft/test-suite/tree/master/MultiSource/Benchmarks/ASC_Sequoia)

### 3.1 Original Version (`IRSmk.ori/`)

* Includes memory tracking extensions

### 3.2 AoS-Optimized Version (`IRSmk.aos/`)

* Uses Array of Structures for layout optimization

**Usage:**

```bash
cd IRSmk

# Build original version
cd IRSmk.ori && make

# Build AoS version
cd ../IRSmk.aos && make

# Run both
../bin/irsmk.ori Input/irsmk_input
../bin/irsmk.aos Input/irsmk_input
```

---

## 4. PolyBench Correlation

**Source:** [PolyBenchC 4.2.1](https://github.com/MatthiasJReisinger/PolyBenchC-4.2.1.git)

* `correlation.c`: Original version with memory tracking
* `correlation_block.c`: Blocked version to enhance data locality

**Usage:**

```bash
cd PolyBenchC-4.2.1
make   # Builds original and blocked versions for all problem sizes

# Run examples
bin/correlation.large
bin/correlation_block.large
bin/correlation.extra
bin/correlation_block.extra
```

---

For questions, please refer to the individual `Makefile` under each subdirectory or contact the project maintainer.

