# mdspass3
`mdspass3`, an extension of mdspass2 with FLTK toolkit,
is an in-house molecular dynamics (MD) visualization and analysis tool developed in Y. Umeno Laboratory
(Institute of Industrial Science, The University of Tokyo).  
It supports several empirical interatomic potentials and provides an interactive GUI for structure editing, MD simulation, visualization, NEB analysis, phonon analysis, and more.

This repository contains the source code, and a pre-built binary for Windows.

---

## Features

- Interactive MD visualization (GUI)
- Support for various empirical potentials  
  (Morse, EAM, ADP, Tersoff, AIREBO (torsional term excluded), etc.)
- NEB (Nudged Elastic Band) analysis
- Phonon dispersion calculation from small unit cells
- Structure editing (add/delete/move atoms, change species)
- Stress analysis (global & local)
- Periodic boundary condition handling
- Configuration measurement tools (distances, bond angles)
- Energy/force-only update mode (no MD)

---

## Directory Structure

```

mdspass3/
├── CONFIG/          # Example configuration file
├── SETDAT/          # Example setting file
├── src/             # Source code, auxiliary libraries & build scripts for macOS
├── CMakeLists.txt   # CMake file
├── forWin/          # Pre-compiled binary for Windows
├── pot/             # Potential files
├── Documents/       # Installation guide and other information
└── README.md        # This file

```

---

## Download and build
**NB: See License below before downloading.**

### Windows (pre-compiled binary)
- https://github.com/yoshiumeno/mdspass3/tree/main/forWin

### Linux/macOS (compile from source)
- Source code: 
 https://github.com/yoshiumeno/mdspass3/tree/main/src

- CMake file: 
 https://github.com/yoshiumeno/mdspass3/tree/main/CMakeLists.txt

### How to compile (Linux/macOS)

##### on Ubuntu:
```  
sudo apt install cmake libfltk1.3-dev libglu1-mesa-dev liblapack-dev libblas-dev libpng-dev libpng-dev  
cmake -S . -B build  
cd build  
make  
```

##### on macOS:
```
brew install cmake fltk mesa-glu lapack openblas libpng
cmake -S . -B build
cd build
make
```

### How to run
Go down to `build/` (Linux or macOS) or `forWin` (Windows) directory and make sure that you have
`CONFIG`, `SETDAT` and `pot/` in the directory.
If not, copy them to the directory.
Simply run `./mdspass3` (Linux or macOS) or `mdspass3.exe` (Windows).


---

## License

At this moment, the use of mdspass3 is limited to our collaborative research projects. If you wish to download or use the software, please make sure to contact Umeno Laboratory in advance.

---

## Contact
Prof. Yoshitaka Umeno
Umeno Laboratory
Institute of Industrial Science, The University of Tokyo

* Website: https://www.cmsm.iis.u-tokyo.ac.jp/index.html
* Email: umeno@iis.u-tokyo.ac.jp 


---