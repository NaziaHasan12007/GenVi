# GenVi: Genomic Eco-Efficiency Analyzer

GenVi is a high-performance genomic analysis pipeline engineered in C. It bridges the gap between complex bioinformatics and environmental sustainability. By integrating real-time IoT sensor data with a multi-threaded alignment engine, GenVi calculates the carbon footprint of genomic sequencing and provides actionable insights to reduce the environmental impact of computational biology.

---

## Key Features

- **Parallel Mapping Engine**  
  High-speed DNA read alignment using pthreads and thread-safe fine-grained locking for hash table indexing.

- **Eco-Efficiency Scoring**  
  A custom algorithm that correlates CPU throughput with hardware temperature and grid carbon intensity.

- **IoT Sensor Integration**  
  Parses real-time energy and temperature data to calculate the true cost of sequencing.

- **Advanced Genomic Metrics**  
  Automated calculation of N50, median coverage, and GC-content using optimized merge sort implementations.

- **Recommendation Engine**  
  Provides intelligent feedback to optimize future sequencing runs for a lower carbon footprint.

---

## Input Specifications

| File Type          | Extension      | Description |
|-------------------|---------------|-------------|
| Reference Genome  | .fna / .fasta | The established sequence used as the mapping backbone |
| Sample Reads      | .fastq        | Raw DNA fragments from the sequencing machine |
| Carbon Intensity  | .csv          | Time-stamped data providing grid carbon footprint (gCO2/kWh) |
| IoT Sensor Data   | .log / .csv   | Real-time monitoring of temperature (°C) and power consumption (Watts) |

---

## Technical Deep-Dive

### Multi-Threaded Alignment

GenVi utilizes a **Seed-and-Extend** strategy. The reference genome is decomposed into k-mers stored in a custom hash table with 4,999,999 buckets.

To maximize performance on multi-core systems:

- **Fine-Grained Locking:** Millions of individual mutexes ensure minimal thread contention during indexing  
- **Thread-Local Buffers:** Mapping results are stored in private memory arrays before being merged into the global coverage map, avoiding performance bottlenecks  

---

### Sustainability Logic

The program calculates the carbon footprint (C) using the following variables:

- **Energy Consumption (E):** Derived from IoT sensor data  
- **PUE (Power Usage Effectiveness):** Adjusted based on machine temperature logs  
- **Grid Intensity (I):** Matched via timestamps from the emissions dataset  

---

## Getting Started

### Prerequisites

- **Compiler:** GCC or Clang (C11 or later)  
- **Library:** pthread (standard on Linux/macOS)  
- **Operating System:** Linux, macOS (M3 optimized), or Windows (via WSL)  

---

### Installation & Execution

```bash
git clone https://github.com/yourusername/genvi.git
cd genvi
gcc -O3 *.c -o genvi -lpthread
./genvi reference.fna sample.fastq emissions.csv sensors.csv
```

---

## Future Enhancements

- GPU Acceleration: Transition Smith-Waterman kernels to CUDA for ~10x better energy efficiency  
- Real-time API Integration: Fetch live grid intensity data from Electricity Maps  
- 2-Bit Encoding: Reduce RAM usage by up to 75% through bit-packing  

---

## Contributors

**Nazia Hassan**  
Software Engineering Student  
University of Dhaka

Copyright (c) 2026 Nazia Hassan

This project is submitted for academic evaluation only.
All rights reserved.
This code may not be reused, redistributed, or modified without permission from the author, except for academic review purposes by instructors and evaluators.
