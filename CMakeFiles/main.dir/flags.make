# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

# compile CXX with /usr/bin/mpicxx
CXX_FLAGS = -I/home/maryhallow/.linuxbrew/Cellar/dealii/HEAD/include -I/home/maryhallow/.linuxbrew/Cellar/dealii/HEAD/include/deal.II/bundled -I/usr/lib/openmpi/include -I/usr/lib/openmpi/include/openmpi -I/home/maryhallow/.linuxbrew/lib/cmake/Pike/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/PikeImplicit/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/PikeBlackBox/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/TrilinosCouplings/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Mesquite/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Didasko/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/CTrilinos/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/ROL/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/MOOCHO/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Rythmos/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/MueLu/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Moertel/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/NOX/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Phalanx/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Intrepid/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Teko/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Stratimikos/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/ShyLU/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/ShyLUCore/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Ifpack2/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Zoltan2/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Anasazi/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Belos/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/ML/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Komplex/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Ifpack/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Pamgen/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Amesos/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Galeri/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/AztecOO/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Pliris/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Isorropia/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/OptiPack/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Thyra/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/ThyraTpetraAdapters/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/ThyraEpetraExtAdapters/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/ThyraEpetraAdapters/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/ThyraCore/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Xpetra/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/EpetraExt/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Tpetra/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/TpetraCore/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/TpetraTSQR/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/TpetraKernels/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/TpetraClassic/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Triutils/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/GlobiPack/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Shards/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Zoltan/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Epetra/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Sacado/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/RTOp/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Teuchos/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/TeuchosKokkosComm/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/TeuchosKokkosCompat/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/TeuchosRemainder/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/TeuchosNumerics/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/TeuchosComm/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/TeuchosParameterList/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/TeuchosCore/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Kokkos/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/KokkosAlgorithms/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/KokkosContainers/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/KokkosCore/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/ThreadPool/../../../include -I/home/maryhallow/.linuxbrew/lib/cmake/Gtest/../../../include -I/home/maryhallow/.linuxbrew/include -I/home/maryhallow/.linuxbrew/opt/superlu_dist/include/superlu_dist -I/home/maryhallow/.linuxbrew/opt/parmetis/include -I/home/maryhallow/.linuxbrew/opt/metis/include    -pedantic -fpic -Wall -Wextra -Wpointer-arith -Wwrite-strings -Wsynth -Wsign-compare -Wswitch -Woverloaded-virtual -Wno-long-long  -Wno-literal-suffix -std=c++11 -Os -w -pipe -march=core2 -Wno-parentheses -Wno-unused-local-typedefs -Og -ggdb -Wa,--compress-debug-sections

CXX_DEFINES = -DDEBUG -DTBB_DO_ASSERT=1 -DTBB_USE_DEBUG

