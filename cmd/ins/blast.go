// Copyright ©2020 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"bufio"
	"bytes"
	"encoding/json"
	"encoding/xml"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"

	"modernc.org/kv"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"

	"github.com/kortschak/ins/blast"
	"github.com/kortschak/ins/internal/store"
)

const (
	xmlFmt = 5
	tabFmt = 6
)

// runBlastTabular runs a BLAST search of the sequences in libs against a database
// constructed from the sequences in query with details from g. The BLAST parameters
// are provided by search. The strings mflags and bflags are passed to makeblastdb
// and blastn as flags without interpretation or checking. If logger is not nil,
// output from the blast executable is written to it.
func runBlastTabular(search blast.Nucleic, query *os.File, libs []library, mx map[string]fragment, mflags, bflags string, logger io.Writer) (*kv.DB, error) {
	search.OutFormat = tabFmt

	opts := &kv.Options{Compare: store.GroupByQueryOrderSubjectLeft}
	hits, err := kv.Create(filepath.Join(filepath.Dir(query.Name()), "forward.db"), opts)
	if err != nil {
		return nil, err
	}

	for _, lib := range libs {
		working, err := workingFile(query, "-working")
		if err != nil {
			return nil, err
		}
		for n := 0; n < maxIters; n++ {
			mkdb, err := blast.MakeDB{DBType: "nucl", In: working, Out: working, ExtraFlags: mflags}.BuildCommand()
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
			search.ExtraFlags = bflags
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

			err = mask(working, lastHits, 'N')
			if err != nil {
				return nil, err
			}

			log.Print("remapping coordinates")
			remapCoords(lastHits, mx)
			const batch = 100
			for i, h := range lastHits {
				if i%batch == 0 {
					err = hits.BeginTransaction()
					if err != nil {
						return nil, err
					}
				}
				key := store.MarshalBlastRecordKey(h)
				// Keep a record of the actual hit purely for
				// correctness auditing; the key has enough
				// information for what we need.
				value, err := json.Marshal(h)
				if err != nil {
					return nil, err
				}
				err = hits.Set(key, value)
				if err != nil {
					return nil, err
				}
				if i%batch == batch-1 || i == len(lastHits)-1 {
					err = hits.Commit()
					if err != nil {
						return nil, err
					}
				}
			}

			err = lib.reset()
			if err != nil {
				return nil, err
			}
		}
	}
	return hits, nil
}

func workingFile(src *os.File, suffix string) (name string, err error) {
	dst, err := os.Create(src.Name() + suffix)
	if err != nil {
		return "", err
	}
	_, err = src.Seek(0, io.SeekStart)
	if err != nil {
		return "", err
	}
	_, err = io.Copy(dst, src)
	if err != nil {
		return "", err
	}
	err = dst.Close()
	if err != nil {
		return "", err
	}
	return dst.Name(), nil
}

// runBlastXML runs a BLAST search of the sequences in libs against a database
// constructed from the sequences in query with details from g. The BLAST parameters
// are provided by search. The strings mflags and bflags are passed to makeblastdb
// and blastn as flags without interpretation or checking. Work is done in workdir
// and if logger is not nil, output from the blast executable is written to it.
func runBlastXML(search blast.Nucleic, g store.BlastRecordKey, query io.Reader, libs []library, workdir, mflags, bflags string, logger io.Writer) ([]*blast.Output, error) {
	search.OutFormat = xmlFmt

	working := filepath.Join(workdir, g.QueryAccVer+"-working")
	mkdb, err := blast.MakeDB{DBType: "nucl", In: "-", Title: g.QueryAccVer, Out: working, ExtraFlags: mflags}.BuildCommand()
	if err != nil {
		return nil, err
	}
	log.Printf("%v < <%s %+d matches>", mkdb, g.QueryAccVer, g.Strand)
	mkdb.Stdin = query
	mkdb.Stdout = logger
	mkdb.Stderr = logger
	err = mkdb.Run()
	if err != nil {
		return nil, err
	}

	var results []*blast.Output
	for _, lib := range libs {
		search.Database = working
		search.Query = lib.name()
		search.ExtraFlags = bflags
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
		var o blast.Output
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
					if strand != g.Strand {
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

var hitID int64

func nextID() int64 {
	hitID++
	return hitID
}

// reportBlast converts BLAST results into blast.Records based on the
// coordinates of a genome region g.
func reportBlast(results []*blast.Output, g store.BlastRecordKey, verbose bool) []blast.Record {
	var remapped []blast.Record
	for _, o := range results {
		for _, it := range o.Iterations {
			if it.QueryId == nil {
				log.Printf("missing query id skipping: %s x %s", g.SubjectAccVer, g.QueryAccVer)
				break
			}

			for _, hit := range it.Hits {
				def := hit.Def
				i := strings.Index(def, " ")
				if i >= 0 {
					def = def[:i]
				}
				desc := strings.Fields(hit.Def[i+1:])
				left, err := strconv.Atoi(desc[0])
				if err != nil {
					panic("invalid left range:" + hit.Def)
				}
				right, err := strconv.Atoi(desc[1])
				if err != nil {
					panic("invalid right range:" + hit.Def)
				}

				if *it.QueryId != g.QueryAccVer {
					log.Printf("dropping remainder of iteration hits: unexpected name: %q want: %q range:%d-%d", *it.QueryId, g.QueryAccVer, left, right)
					for i, hsp := range hit.Hsps {
						log.Printf("\t%q:%d-%d %q:%d-%d",
							hit.Def, hsp.HitFrom, hsp.HitTo,
							*it.QueryDef, hsp.QueryFrom, hsp.QueryTo,
						)
						if !verbose {
							n := len(hit.Hsps[i:])
							if n > 2 {
								log.Printf("\tnot logging %d additional dropped hits", n-1)
								break
							}
						}
					}
					break
				}

				// TODO: Handle strand mismatches in a similar way to
				// how we handle ID mismatches. This is slightly more
				// complicated since we don't know strand until we get
				// into the HSPs.

				uid := nextID()
				for _, hsp := range hit.Hsps {
					// Convert to 0-based indexing.
					hsp.QueryFrom--
					hsp.HitFrom--

					// Remap coordinates onto original subject.
					hsp.HitFrom += left
					hsp.HitTo += left

					remapped = append(remapped, blast.Record{
						QueryAccVer: g.QueryAccVer,
						QueryStart:  hsp.QueryFrom,
						QueryEnd:    hsp.QueryTo,

						SubjectAccVer: g.SubjectAccVer,
						SubjectStart:  hsp.HitFrom,
						SubjectEnd:    hsp.HitTo,

						Strand: g.Strand,

						PctIdentity:     100 * float64(*hsp.HspIdentity) / float64(*hsp.AlignLen),
						AlignmentLength: *hsp.AlignLen,
						Mismatches:      *hsp.AlignLen - *hsp.HspIdentity,
						GapOpens:        *hsp.HspGaps,
						EValue:          hsp.EValue,
						BitScore:        hsp.BitScore,

						UID: uid,
					})
				}
			}
		}
	}

	return remapped
}

// mask writes a masked copy of the genome in the src file based on the given
// blast hits. Regions that are masked are replaced with the masked alphabet.Letter.
func mask(path string, hits []blast.Record, masked alphabet.Letter) error {
	log.Printf("masking %s", path)
	src, err := os.Open(path)
	if err != nil {
		return err
	}
	defer src.Close()

	dst, err := ioutil.TempFile(filepath.Dir(path), filepath.Base(path)+"-*")
	if err != nil {
		return err
	}
	defer dst.Close()

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
		fmt.Fprintf(dst, "%60a\n", seq)
	}
	err = sc.Error()
	if err != nil {
		return err
	}
	err = dst.Sync()
	if err != nil {
		return err
	}
	return os.Rename(dst.Name(), path)
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
