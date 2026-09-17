# compile
    cmake -S ./
    make

We get bin/mce and bin/kcc now.

# data format and command

## each edge per line

    bin/mce noUVM -f_txt dataFile
    bin/mce noUVM -f_txt dataFile --print-cliques
    bin/mce noUVM -f_txt dataFile --output-cliques cliques.txt
    bin/kcc noUVM -f_txt dataFile -k k
    bin/KCCparallel noUVM -f_txt dataFile -k k -steal

## The first line is the number of nodes and edges, then each edge per line

    bin/mce -f_txt dataFile
    bin/kcc -f_txt dataFile -k k
    bin/MCEparallel -f_txt dataFile -k k -t 10 -steal

# local experiment configuration

The sequential `bin/mce` path used by `../baselines/run.sh` counts maximal
cliques of size at least 3. Every counted clique is constructed explicitly and
retained in memory as a sorted list of original input vertex labels; there is
no count-only or closed-form counting path. The executable prints
`stored_cliques` and requires it to equal `Mclique`, plus
`clique_storage:all` to make this contract visible.

`--print-cliques` prints the already-retained cliques after the timed
enumeration. `--output-cliques FILE` saves the same retained list with one
space-separated clique per line and no header. Either, both, or neither option
may be used; storage always occurs. Retaining every identity can consume
substantial memory on graphs with millions of maximal cliques.

`../../tests/run_external_baselines_smoke.sh` includes an exhaustive
differential check against brute force for all 33,868 labeled simple graphs
with at most six vertices, in addition to saved/printed-list fixtures.

The local file-size helper is also safe in Release/NDEBUG builds; upstream's
assert-only `stat` call was not.
