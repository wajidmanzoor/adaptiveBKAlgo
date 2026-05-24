# Maximal Clique Enumeration with Hybrid Branching and Early Termination

This repository contains the codes and datasets for maximal clique enumeration employed with the edge-oriented branching strategy and the early-termination technique. 

## Setup

To run the codes correctly, please ensure that (1) you have installed `cmake (minimum version 3.16)` on your computer or server since we need to use cmake to generate `Makefile` to compile our codes and (2) the graph data should contain no self-loops and no duplicate edges. 

## Step 0 - Download and compile the codes

To download the codes and datasets, there are two options: (1) you can download HBBMC.zip and unzip it or (2) you can download each file in `src/` and `dataset/`. Remember to keep the forlder's structure, shown as follows. 

```bash
.(HBBMC)
├── dataset
│   ├── dblp.clean
│   └── youtube.clean
└── src
    ├── CMakeLists.txt
    ├── edge_oriented.cpp
    ├── edge_oriented.h
    ├── enum.cpp
    ├── enum.h
    ├── graph.hpp
    ├── hash.hpp
    ├── tool.cpp
    └── tool.h
``` 

Run the following commands to compile our codes. 

```bash
cd src/
mkdir build && cd build/
cmake ..
make
```

After running the above commands, there will be an executable file called `MCE`. 

## Step 1 - Enumeration Procedure

To enumerate all maximal cliques (MCs) in the graph, the running command is: 

```bash
./MCE /PATH_TO_DATA t
```

For example, 

```bash
./MCE ../../dataset/dblp.clean 3     # Enumerate MCs in `dblp` with early-termination in 3-plex
./MCE ../../dataset/youtube.clean 2  # Enumerate MCs in `youtube` with early-termination in 2-plex
```

For `dblp` dataset, it outputs

```
Reading the graph.
Finish reading the graph. |V| = 317080, |E| = 1049866
#mc = 257551, runtime = 163.666 ms
```

For `youtube` dataset, it outputs

```
Reading the graph.
Finish reading the graph. |V| = 1134890, |E| = 2987624
#mc = 3265956, runtime = 1459.230 ms
```

