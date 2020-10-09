// Copyright ©2020 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"bufio"
	"bytes"
	"encoding/xml"
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	xmlblast "github.com/biogo/ncbi/blast"

	"github.com/kortschak/ins/blast"
)

const (
	xmlFmt = 5
	tabFmt = 6
)

// runBlastTabular runs a BLAST search of the sequences in libs against a database
// constructed from the sequences in query with details from g. The BLAST parameters
// are provided by search. If logger is not nil, output from the blast executable is
// written to it.
func runBlastTabular(search blast.Nucleic, query *os.File, libs []library, logger io.Writer) ([]blast.Record, error) {
	search.OutFormat = tabFmt

	var hits []blast.Record
	for _, lib := range libs {
		for n := 0; n < maxIters; n++ {
			working, err := mask(query, query.Name()+"-working", hits, 'N')
			if err != nil {
				return nil, err
			}

			mkdb, err := blast.MakeDB{DBType: "nucl", In: working, Out: working}.BuildCommand()
			if err != nil {
				return nil, err
			}
			log.Print(mkdb)
			mkdb.Stdout = logger
			mkdb.Stderr = logger
			err = mkdb.Run()
			if err != nil {
				return nil, err
			}

			search.Database = working
			search.Query = lib.name()
			blastn, err := search.BuildCommand()
			if err != nil {
				return nil, err
			}

			log.Print(blastn)
			blastn.Stdin = lib.stream()
			blastn.Stderr = logger
			stdout, err := blastn.StdoutPipe()
			if err != nil {
				return nil, err
			}
			err = blastn.Start()
			if err != nil {
				return nil, err
			}

			lastHits, err := blast.ParseTabular(stdout, n)
			if err != nil {
				return nil, err
			}
			log.Printf("blast iteration %d found %d new matches", n, len(lastHits))

			err = blastn.Wait()
			if err != nil {
				return nil, err
			}

			if len(lastHits) == 0 {
				break
			}
			hits = append(hits, lastHits...)

			err = lib.reset()
			if err != nil {
				return nil, err
			}
		}
	}
	return hits, nil
}

// runBlastXML runs a BLAST search of the sequences in libs against a database
// constructed from the sequences in query with details from g. The BLAST parameters
// are provided by search. Work is done in workdir and if logger is not nil, output
// from the blast executable is written to it.
func runBlastXML(search blast.Nucleic, g blastRecordGroup, query io.Reader, libs []library, workdir string, logger io.Writer) ([]*xmlblast.Output, error) {
	search.OutFormat = xmlFmt

	working := filepath.Join(workdir, g.QueryAccVer+"-working")
	mkdb, err := blast.MakeDB{DBType: "nucl", In: "-", Title: g.QueryAccVer, Out: working}.BuildCommand()
	if err != nil {
		return nil, err
	}
	log.Printf("%v < <%s %+d matches>", mkdb, g.QueryAccVer, g.strand)
	mkdb.Stdin = query
	mkdb.Stdout = logger
	mkdb.Stderr = logger
	err = mkdb.Run()
	if err != nil {
		return nil, err
	}

	var results []*xmlblast.Output
	for _, lib := range libs {
		search.Database = working
		search.Query = lib.name()
		blastn, err := search.BuildCommand()
		if err != nil {
			return nil, err
		}

		log.Print(blastn)
		blastn.Stdin = lib.stream()
		blastn.Stderr = logger
		stdout, err := blastn.StdoutPipe()
		if err != nil {
			return nil, err
		}
		err = blastn.Start()
		if err != nil {
			return nil, err
		}

		dec := xml.NewDecoder(stdout)
		var o xmlblast.Output
		err = dec.Decode(&o)
		if err != nil {
			return nil, err
		}

		err = blastn.Wait()
		if err != nil {
			return nil, err
		}

		i := 0
		for _, it := range o.Iterations {
			if len(it.Hits) == 0 || it.QueryId == nil || *it.QueryId != g.QueryAccVer {
				continue
			}

			for j, hit := range it.Hits {
				k := 0
				for _, hsp := range hit.Hsps {
					qStrand := 1
					if hsp.QueryTo < hsp.QueryFrom {
						qStrand = -1
					}
					hStrand := 1
					if hsp.HitTo < hsp.HitFrom {
						hStrand = -1
					}
					strand := int8(qStrand * hStrand)
					if strand != g.strand {
						continue
					}
					hit.Hsps[k] = hsp
					k++
				}
				it.Hits[j].Hsps = hit.Hsps[:k]
			}

			o.Iterations[i] = it
			i++
		}
		o.Iterations = o.Iterations[:i]

		results = append(results, &o)
	}
	return results, nil
}

// reportBlast converts BLAST results into blast.Records based on the
// coordinates of a genome region g.
func reportBlast(results []*xmlblast.Output, g blastRecordGroup) []blast.Record {
	var remapped []blast.Record
	for _, o := range results {
		for _, it := range o.Iterations {
			for _, hit := range it.Hits {
				id := hit.Def
				i := strings.IndexAny(id, " \t")
				if i >= 0 {
					id = id[:i]
				}
				id = strings.TrimSuffix(id, fmt.Sprintf("_%d_%d", g.left, g.right))
				if id != g.SubjectAccVer {
					panic("unexpected name: " + id)
				}

				for _, hsp := range hit.Hsps {
					// Convert to 0-based indexing.
					hsp.QueryFrom--
					hsp.HitFrom--

					// Remap coordinates onto original subject.
					hsp.HitFrom += g.left
					hsp.HitTo += g.left

					remapped = append(remapped, blast.Record{
						QueryAccVer: g.QueryAccVer,
						QueryStart:  hsp.QueryFrom,
						QueryEnd:    hsp.QueryTo,

						SubjectAccVer: g.SubjectAccVer,
						SubjectStart:  hsp.HitFrom,
						SubjectEnd:    hsp.HitTo,

						Strand: g.strand,

						PctIdentity:     100 * float64(*hsp.HspIdentity) / float64(*hsp.AlignLen),
						AlignmentLength: *hsp.AlignLen,
						Mismatches:      *hsp.AlignLen - *hsp.HspIdentity,
						GapOpens:        *hsp.HspGaps,
						EValue:          hsp.EValue,
						BitScore:        hsp.BitScore,
					})
				}
			}
		}
	}

	return remapped
}

// mask writes a masked copy of the genome in the src file based on the given
// blast hits. Regions that are masked are replaced with the masked alphabet.Letter.
func mask(src *os.File, dstPath string, hits []blast.Record, masked alphabet.Letter) (name string, err error) {
	log.Printf("masking %s", src.Name())
	_, err = src.Seek(0, io.SeekStart)
	if err != nil {
		return "", err
	}

	working, err := os.Create(dstPath)
	if err != nil {
		return "", err
	}
	defer working.Close()

	if len(hits) == 0 {
		_, err = io.Copy(working, src)
		if err != nil {
			return working.Name(), err
		}
		return working.Name(), working.Sync()
	}

	hitsOf := make(map[string][]blast.Record)
	for _, h := range hits {
		hitsOf[h.SubjectAccVer] = append(hitsOf[h.SubjectAccVer], h)
	}

	sc := seqio.NewScanner(fasta.NewReader(src, linear.NewSeq("", nil, alphabet.DNAredundant)))
	for sc.Next() {
		seq := sc.Seq().(*linear.Seq)
		for _, h := range hitsOf[seq.ID] {
			// Blast reports minus strand matches by inverting the coordinates.
			if h.SubjectEnd < h.SubjectStart {
				h.SubjectStart, h.SubjectEnd = h.SubjectEnd, h.SubjectStart
			}
			for i := h.SubjectStart; i < h.SubjectEnd; i++ {
				seq.Seq[i-seq.Offset] = masked
			}
		}
		fmt.Fprintf(working, "%60a\n", seq)
	}
	err = sc.Error()
	if err != nil {
		return working.Name(), err
	}

	return working.Name(), working.Sync()
}

// detail is the class and length of a repeat type.
type detail struct {
	class  string
	length int
}

// libDetails returns the details of the repeats in lib.
func libDetails(lib []library) (map[string]detail, error) {
	details := make(map[string]detail)
	for _, l := range lib {
		var r io.Reader
		switch l := l.(type) {
		case *stream:
			err := l.reset()
			if err != nil {
				return nil, err
			}
			r = l.stream()
		case filename:
			f, err := os.Open(l.name())
			if err != nil {
				return nil, err
			}
			defer f.Close()
			r = f
		default:
			panic("unknown library type")
		}
		sc := bufio.NewScanner(r)
		sc.Split(func(data []byte, atEOF bool) (advance int, token []byte, err error) {
			if atEOF && len(data) == 0 {
				return 0, nil, nil
			}
			if i := bytes.IndexByte(data, '\n'); i >= 0 {
				return i + 1, data[:i+1], nil
			}
			if atEOF {
				return len(data), data, nil
			}
			return 0, nil, nil
		})

		var (
			name   string
			rec    detail
			offset int64
		)
		for sc.Scan() {
			b := bytes.TrimSpace(sc.Bytes())
			if len(b) == 0 {
				continue
			}
			if b[0] == '>' {
				if name != "" {
					details[name] = rec
					rec = detail{}
				}
				lenID := bytes.IndexAny(b, " \t")
				if lenID < 0 {
					name = string(b[1:])
					rec.class = ""
				} else {
					name = string(b[1:lenID])
					rec.class = string(bytes.Fields(b[lenID+1:])[0])
				}
				if _, exists := details[name]; exists {
					return nil, fmt.Errorf("duplicate sequence identifier %s at %d", name, offset)
				}
			} else {
				rec.length += len(b)
			}
			offset += int64(len(sc.Bytes()))
		}
		if name != "" {
			details[name] = rec
		}
		err := sc.Err()
		if err != nil {
			return nil, err
		}
	}
	return details, nil
}

type library interface {
	name() string
	stream() io.Reader
	reset() error
}

func newStream(s []string) ([]library, error) {
	names := uniq(s)
	st := &stream{files: make([]*os.File, len(names))}
	for i, n := range names {
		f, err := os.Open(n)
		if err != nil {
			return nil, err
		}
		st.files[i] = f
	}
	err := st.reset()
	if err != nil {
		return nil, err
	}
	return []library{st}, nil
}

type stream struct {
	multi io.Reader

	files []*os.File
}

func (s *stream) name() string      { return "-" }
func (s *stream) stream() io.Reader { return s.multi }
func (s *stream) reset() error {
	s.multi = nil
	r := make([]io.Reader, len(s.files))
	for i, f := range s.files {
		_, err := f.Seek(0, io.SeekStart)
		if err != nil {
			return err
		}
		r[i] = f
	}
	s.multi = io.MultiReader(r...)
	return nil
}

func filenames(s []string) []library {
	f := make([]library, len(s))
	for i, v := range s {
		f[i] = filename(v)
	}
	return f
}

type filename string

func (f filename) name() string      { return string(f) }
func (f filename) reset() error      { return nil }
func (f filename) stream() io.Reader { return nil }

func uniq(s []string) []string {
	sort.Strings(s)
	i := 0
	for _, v := range s {
		if v != s[i] {
			i++
			s[i] = v
		}
	}
	return s[:i+1]
}

type interval struct {
	parent     string
	start, end int
}
