// Copyright ©2020 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"fmt"
	"io"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"

	"github.com/kortschak/ins/blast"
)

// split splits the fasta sequence read from src into fragments that are no longer
// than max but segmenting into fragments that are goal long. It writes the coordinates
// of the sequence relative to the original in the first three space separated fields
// of the fasta description and returns a map containing a look-up table from the
// generated sequences to the parent and coordinates.
func split(dst io.Writer, src io.Reader, goal, max int) (map[string]fragment, error) {
	frags := make(map[string]fragment)
	sc := seqio.NewScanner(fasta.NewReader(src, linear.NewSeq("", nil, alphabet.DNA)))
	i := 1
	for sc.Next() {
		pos := 0
		seq := sc.Seq().(*linear.Seq)
		id := seq.ID
		desc := seq.Desc
		for seq.Len() > max {
			tmp := *seq
			n := min(len(tmp.Seq), goal)
			tmp.Seq = tmp.Seq[:n]
			tmp.ID = fmt.Sprintf("%s_%d", id, i)
			tmp.Desc = fmt.Sprintf("%s %d %d %s", id, pos, pos+n, desc)
			if _, ok := frags[tmp.ID]; ok {
				return nil, fmt.Errorf("non-unique sequence id in input: %q", id)
			}
			frags[tmp.ID] = fragment{parent: id, start: pos, end: pos + n}
			fmt.Fprintf(dst, "%60a\n", &tmp)
			seq.Seq = seq.Seq[n:]
			pos += n
			i++
		}
		seq.ID = fmt.Sprintf("%s_%d", id, i)
		seq.Desc = fmt.Sprintf("%s %d %d %s", id, pos, pos+seq.Len(), desc)
		if _, ok := frags[seq.ID]; ok {
			return nil, fmt.Errorf("non-unique sequence id in input: %q", id)
		}
		frags[seq.ID] = fragment{parent: id, start: pos, end: pos + seq.Len()}
		fmt.Fprintf(dst, "%60a\n", seq)
	}
	if err := sc.Error(); err != nil {
		return nil, fmt.Errorf("error during sequence read: %w", err)
	}
	return frags, nil
}

// remapCoords adjusts hits so that subjects (genome sequence) are mapped against
// the original un-fragmented genome sequence consumed by split. It then sorts
// hits by strand, repeat type, position and BLAST bitscore.
func remapCoords(hits []blast.Record, frags map[string]fragment) {
	for i, r := range hits {
		iv := frags[r.SubjectAccVer]
		r.SubjectAccVer = iv.parent
		r.SubjectStart += iv.start
		r.SubjectEnd += iv.start
		hits[i] = r
	}
}

type fragment struct {
	parent     string
	start, end int
}

// merge takes a sorted set of hits and groups them into individual groups based
// on proximity. If adjacent hits are within near, they are grouped.
func merge(hits []blast.Record, near int) (groups []blastRecordGroup) {
	if len(hits) == 0 {
		return nil
	}

	r := hits[0]
	groups = []blastRecordGroup{{
		QueryAccVer:   r.QueryAccVer,
		SubjectAccVer: r.SubjectAccVer,
		left:          min(r.SubjectStart, r.SubjectEnd),
		right:         max(r.SubjectStart, r.SubjectEnd),
		strand:        r.Strand,
		recs:          []blast.Record{r},
	}}
	for _, r := range hits[1:] {
		left := min(r.SubjectStart, r.SubjectEnd)
		right := max(r.SubjectStart, r.SubjectEnd)

		last := &groups[len(groups)-1]
		if left-last.right <= near && r.Strand == last.strand && r.SubjectAccVer == last.SubjectAccVer && r.QueryAccVer == last.QueryAccVer {
			last.recs = append(last.recs, r)
			last.right = max(last.right, right)
			continue
		}

		groups = append(groups, blastRecordGroup{
			QueryAccVer:   r.QueryAccVer,
			SubjectAccVer: r.SubjectAccVer,
			left:          left,
			right:         right,
			strand:        r.Strand,
			recs:          []blast.Record{r},
		})
	}

	return groups
}

// blastREcord group is a set of blast hits that are grouped by proximity.
type blastRecordGroup struct {
	QueryAccVer   string
	SubjectAccVer string
	left          int
	right         int
	strand        int8

	recs []blast.Record
}

// bySubjectLeft satisfies the sort.Interface, ordering by strand, query name,
// subject name, subject position and BLAST bitscore.
type bySubjectLeft []blast.Record

func (r bySubjectLeft) Len() int {
	return len(r)
}
func (r bySubjectLeft) Less(i, j int) bool {
	// Separate strands, (+) first.
	if r[i].Strand != r[j].Strand {
		return r[i].Strand > r[j].Strand
	}

	// Group elements of the same type.
	if r[i].QueryAccVer != r[j].QueryAccVer {
		return r[i].QueryAccVer < r[j].QueryAccVer
	}

	// Sort by left position, with higher scoring matches first.
	if r[i].SubjectAccVer != r[j].SubjectAccVer {
		return r[i].SubjectAccVer < r[j].SubjectAccVer
	}
	leftI := min(r[i].SubjectStart, r[i].SubjectEnd)
	leftJ := min(r[j].SubjectStart, r[j].SubjectEnd)
	if leftI < leftJ {
		return true
	}
	return leftI == leftJ && r[i].BitScore > r[j].BitScore
}
func (r bySubjectLeft) Swap(i, j int) {
	r[i], r[j] = r[j], r[i]
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}
