# RACES
RACES is an Advanced Cancer Evolution Simulator.

RACES is framework for simulating cancer genomic evolution. It supports:
-   cell strand competition and spacial simulations
-   epigenetic switches
-   copy-number variations
-   SBS signature support

RACES consists in 4 main programs:
-   `build_context_index` creates an index for the mutational contexts (i.e., consecutive triplets of nucleotides) in a genome. In normal conditions, the context index for a genome must be created once for all and the resulting file can be used many times.
-   `species_sim` performs cell spacial simulation by using mutant evolutional information. It supports epigenetic switching.
-   `tissue_sampler` samples the results of `species_sim`.
-   `mutations_sim` builds the phylogenetic forest for the sampled cells, and places mutations (only SNVs and CNAs are supported at the moment) on the tree. It can:
    *   produce simulated reads and generate the correspoding SAM file for them;
    *   collect SNVs and CNAs statistical data and save them in CSV files.

## Required Packages
-   CMake
-   a C++17 compiler
-   Boost::program_options
-   nlohmann/json library

### Optional Packages
-   Boost::unit_test_framework (it enables C++ code testing)
-   SDL_ttf and SDL_image libraries (they enable simulation plotting)
-   Boost::Python (it enables Python bindings)

## Setup

### GNU/Linux and macOS

Once the [required packages](#required-packages) have been installed, use the following commands in a shell

```bash
git clone https://github.com/albertocasagrande/RACES.git
cd RACES
cmake .
make -j
```

### Windows

TBD

## Usage

TBD

## License

Copyright (c) 2023 
Alberto Casagrande <[alberto.casagrande@uniud.it](mailto:alberto.casagrande@uniud.it)>

MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
 
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


