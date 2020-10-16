// Copyright ©2020 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// ins is a genomic repeat identification tool. It finds interspersed repeat
// elements in a genome from a repeat library, makes an N-masked copy of the
// genome and gives a list of found elements either in JSON format or GTF.
package main

import (
	"bufio"
	"bytes"
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"os"
	"path/filepath"
	"runtime"
	"sort"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/featio/gff"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/biogo/seq/linear"
	"github.com/biogo/hts/fai"
	"github.com/biogo/store/interval"

	"github.com/kortschak/ins/blast"
)

var (
	// blastnModes are the first pass BLAST parameters.
	blastnModes = map[string]blast.Nucleic{
		"sensitive": {NumAlignments: 1e7, SearchSpace: 1e6, EValue: 3e-5, Threads: runtime.NumCPU(), Reward: 3, Penalty: -4, GapOpen: 30, GapExtend: 6, XdropUngap: 80, XdropGap: 130, XdropGapFinal: 150, WordSize: 9, ParseDeflines: true, OutFormat: tabFmt},
		"normal":    {NumAlignments: 1e7, SearchSpace: 1e6, EValue: 2e-5, Threads: runtime.NumCPU(), Reward: 3, Penalty: -4, GapOpen: 30, GapExtend: 6, XdropUngap: 80, XdropGap: 130, XdropGapFinal: 150, WordSize: 10, ParseDeflines: true, OutFormat: tabFmt},
		"rough":     {NumAlignments: 1e7, SearchSpace: 1e6, EValue: 1e-5, Threads: runtime.NumCPU(), Reward: 3, Penalty: -4, GapOpen: 30, GapExtend: 6, XdropUngap: 80, XdropGap: 130, XdropGapFinal: 150, WordSize: 11, ParseDeflines: true, OutFormat: tabFmt},
	}

	// realign is the reciprocal hit pass BLAST parameters.
	realign = blast.Nucleic{NumAlignments: 1e7, SearchSpace: 1e6, EValue: 1e-5, Threads: runtime.NumCPU(), Reward: 3, Penalty: -4, GapOpen: 30, GapExtend: 6, XdropUngap: 80, XdropGap: 150, XdropGapFinal: 150, WordSize: 11, ParseDeflines: true, Dust: &blast.Dust{Filter: true}, SoftMask: true, OutFormat: xmlFmt}

	// classify is the repeat classification BLAST parameters.
	//nolint:deadcode,unused,varcheck
	classify = blast.Nucleic{NumAlignments: 1e7, SearchSpace: 1e6, EValue: 1e-5, Threads: runtime.NumCPU(), Reward: 3, Penalty: -4, GapOpen: 30, GapExtend: 6, XdropUngap: 80, XdropGap: 150, XdropGapFinal: 150, WordSize: 11, ParseDeflines: true, Dust: &blast.Dust{Filter: true}, SoftMask: true, OutFormat: tabFmt}
)

const (
	// Maximum number of first pass BLAST searches.
	maxIters = 100
	// Optimal fragment length to split genome into.
	optFragmentLen = 100000
	// Maximum fragment length to split genome into.
	maxFragmentLen = 150000
)

// near is the distance between hits of the same type on the same stand that will
// be merged into one hit.
const near = 30

func main() {
	var libs sliceValue
	in := flag.String("query", "", "specify query sequence file (required)")
	flag.Var(&libs, "lib", "specify the search libraries (required - may be present more than once)")
	mode := flag.String("mode", "normal", "specify search mode")
	jsonOut := flag.Bool("json", false, "specify json format for feature output")
	cull := flag.Bool("cull", true, "specify to discard lower scoring nested features")
	verbose := flag.Bool("verbose", false, "specify verbose logging")
	pool := flag.Bool("pool", true, "specify to pool all libraries into a single search")
	threads := flag.Int("cores", 0, "specify the maximum number of cores for blast searches (<=0 is use all cores)")
	work := flag.Bool("work", false, "specify to keep temporary files")

	flag.Usage = func() {
		fmt.Fprintf(flag.CommandLine.Output(), `Usage of %[1]s:
  $ %[1]s [options] -lib <library.fa> [-lib <library.fa> ...] -query <seq.fa> >out.gtf 2>out.log

Options:
`, os.Args[0])
		flag.PrintDefaults()
	}

	flag.Parse()

	if *in == "" || len(libs) == 0 {
		flag.Usage()
		os.Exit(2)
	}

	search, ok := blastnModes[*mode]
	if !ok {
		log.Fatalf("unknown search mode: %q", *mode)
	}
	if *threads > 0 {
		search.Threads = min(*threads, search.Threads)
	}

	log.Println(os.Args)
	var logger io.WriteCloser
	if *verbose {
		logger = logCapture()
		defer logger.Close()
	}

	tmpDir, err := ioutil.TempDir("", "ins-tmp-*")
	if err != nil {
		log.Fatal(err)
	}
	log.Printf("working in %s", tmpDir)
	if *work {
		log.Println("keeping work")
	} else {
		defer func() {
			os.RemoveAll(tmpDir)
		}()
	}

	query, err := os.Open(*in)
	if err != nil {
		log.Fatal(err)
	}
	defer query.Close()

	frags, err := ioutil.TempFile(tmpDir, "dbmx-*")
	if err != nil {
		log.Fatal(err)
	}

	log.Println("indexing query")
	qidx, err := fai.NewIndex(query)
	if err != nil {
		log.Fatal(err)
	}
	_, err = query.Seek(0, io.SeekStart)
	if err != nil {
		log.Fatal(err)
	}

	log.Println("splitting query")
	mx, err := split(frags, query, optFragmentLen, maxFragmentLen)
	if err != nil {
		log.Fatal(err)
	}
	err = frags.Sync()
	if err != nil {
		log.Fatal(err)
	}

	var libraries []library
	libs = uniq(libs)
	if len(libs) > 1 && *pool {
		libraries, err = newStream(libs)
		if err != nil {
			log.Fatal(err)
		}
	} else {
		libraries = filenames(libs)
	}
	hits, err := runBlastTabular(search, frags, libraries, logger)
	if err != nil {
		log.Fatal(err)
	}

	remapCoords(hits, mx)

	var remappedHits []blast.Record
	qfa := fai.NewFile(query, qidx)
	groups := merge(hits, near)
	var buf bytes.Buffer
	for i, g := range groups {
		seq, err := qfa.SeqRange(g.SubjectAccVer, g.left, g.right)
		if err != nil {
			log.Fatal(err)
		}
		b, err := ioutil.ReadAll(seq)
		if err != nil {
			log.Fatal(err)
		}
		s := linear.NewSeq(fmt.Sprintf("%s_%d_%d", g.SubjectAccVer, g.left, g.right), alphabet.BytesToLetters(b), alphabet.DNAredundant)
		s.Desc = fmt.Sprintf("%d %d %s %+d", g.left, g.right, g.QueryAccVer, g.strand)
		fmt.Fprintf(&buf, "%60a\n", s)

		if i == len(groups)-1 || g.QueryAccVer != groups[i+1].QueryAccVer {
			var libraries []library
			if len(libs) > 1 && *pool {
				libraries, err = newStream(libs)
				if err != nil {
					log.Fatal(err)
				}
			} else {
				libraries = filenames(libs)
			}

			hits, err := runBlastXML(realign, g, &buf, libraries, tmpDir, logger)
			if err != nil {
				log.Fatal(err)
			}

			remappedHits = append(remappedHits, reportBlast(hits, g)...)

			buf.Reset()
		}
	}

	if *cull {
		remappedHits = cullContained(remappedHits)
	}

	sort.Sort(bySubjectPosition(remappedHits))

	if *jsonOut {
		enc := json.NewEncoder(os.Stdout)
		for _, r := range remappedHits {
			err = enc.Encode(r)
			if err != nil {
				log.Fatalf("failed to write feature: %v", err)
			}
		}
	} else {
		details, err := libDetails(libraries)
		if err != nil {
			log.Fatalf("failed to get feature lengths: %v", err)
		}
		enc := gff.NewWriter(os.Stdout, 60, true)
		for _, r := range remappedHits {
			if r.Strand < 0 {
				r.SubjectStart, r.SubjectEnd = r.SubjectEnd, r.SubjectStart
			}
			repeat := details[r.QueryAccVer]
			_, err = enc.Write(&gff.Feature{
				SeqName:    r.SubjectAccVer,
				Source:     "ins",
				Feature:    "repeat",
				FeatStart:  r.SubjectStart,
				FeatEnd:    r.SubjectEnd,
				FeatScore:  &r.BitScore,
				FeatStrand: seq.Strand(r.Strand),
				FeatFrame:  gff.NoFrame,
				FeatAttributes: gff.Attributes{{
					Tag:   "Repeat",
					Value: fmt.Sprintf("%s %s %d %d %d", r.QueryAccVer, repeat.class, r.QueryStart+1, r.QueryEnd, repeat.length-r.QueryEnd),
				}},
			})
			if err != nil {
				log.Fatalf("failed to write feature: %v", err)
			}
		}
	}

	masked, err := mask(query, filepath.Base(query.Name())+"-masked.fasta", remappedHits, 'N')
	if err != nil {
		log.Fatal(err)
	}
	log.Printf("masked sequence in %s", masked)
}

// cullContained returns a copy of hits with all hits that are completely contained by
// a higher scoring hit removed.
func cullContained(hits []blast.Record) []blast.Record {
	var tree interval.IntTree
	for i, r := range hits {
		err := tree.Insert(subjectInterval{uid: uintptr(i), Record: r}, true)
		if err != nil {
			panic(fmt.Sprint(err))
		}
	}
	tree.AdjustRanges()
	var culled []blast.Record
outer:
	for _, r := range hits {
		o := tree.Get(subjectInterval{Record: r})
		for _, h := range o {
			if h.(subjectInterval).BitScore > r.BitScore {
				continue outer
			}
		}
		culled = append(culled, r)
	}
	return culled
}

type subjectInterval struct {
	uid uintptr
	blast.Record
}

// Overlap returns whether the b interval completely contains i.
func (i subjectInterval) Overlap(b interval.IntRange) bool {
	left := min(i.SubjectStart, i.SubjectEnd)
	right := max(i.SubjectStart, i.SubjectEnd)
	return b.Start <= left && right <= b.End
}
func (i subjectInterval) ID() uintptr { return i.uid }
func (i subjectInterval) Range() interval.IntRange {
	return interval.IntRange{Start: min(i.SubjectStart, i.SubjectEnd), End: max(i.SubjectStart, i.SubjectEnd)}
}

type bySubjectPosition []blast.Record

func (r bySubjectPosition) Len() int { return len(r) }
func (r bySubjectPosition) Less(i, j int) bool {
	if r[i].SubjectAccVer < r[j].SubjectAccVer {
		return true
	}
	return r[i].SubjectAccVer == r[j].SubjectAccVer && min(r[i].SubjectStart, r[i].SubjectEnd) < min(r[j].SubjectStart, r[j].SubjectEnd)
}
func (r bySubjectPosition) Swap(i, j int) { r[i], r[j] = r[j], r[i] }

// sliceValue is a multi-value flag value.
type sliceValue []string

// Set adds the string to the sliceValue.
func (s *sliceValue) Set(v string) error {
	*s = append(*s, v)
	return nil
}

// String satisfies the flag.Value interface.
func (s *sliceValue) String() string {
	return fmt.Sprintf("%q", []string(*s))
}

// logCapture returns an io.WriteCloser that pipes writes to the default log logger.
func logCapture() io.WriteCloser {
	r, w := io.Pipe()
	go func() {
		sc := bufio.NewScanner(r)
		for sc.Scan() {
			if len(bytes.TrimSpace(sc.Bytes())) == 0 {
				continue
			}
			log.Printf("\t%s", sc.Bytes())
		}
		err := sc.Err()
		if err != nil && err != io.EOF {
			_ = w.CloseWithError(err)
		}
	}()
	return w
}
