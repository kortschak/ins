# `ins`

The `ins` tool is a repeat identification/annotation tool.

`ins` takes a query sequence (usually a genome) and a library or set of FASTA libraries of know repeat sequences. It uses nucleotide BLAST+ to find instances of the repeats in the library within the query sequence and outputs a copy of the sequence with the identified repeats masked with N. The locations of repeats are output on standard output in either GTF format or as a JSON stream. Logging is made to standard error.

## Installation

The `ins` program can be installed after [installing the Go programming language](https://golang.org/doc/install).

```
$ go get github.com/kortschak/cmd/ins
```

`ins` requires that the NCBI+ BLAST distribution is installed on your system and that `blastn` and `makeblastdb` are in your `$PATH`.

## Usage

```
$ ins [options] -lib <library.fa> [-lib <library.fa> ...] -query <seq.fa> >out.gtf 2>out.log
```
or
```
$ ins [options] -json -lib <library.fa> [-lib <library.fa> ...] -query <seq.fa> >out.json 2>out.log
```
