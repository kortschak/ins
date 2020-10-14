// Copyright ©2020 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// cull is a tool to remove lower scoring features from a GFF file.
// It discards features that are completely contained within a higher scoring
// feature. Features without a score are not considered but retained in the
// set of features.
//
// usage: cull < infile.gff > outfile.gff
package main

import (
	"flag"
	"fmt"
	"log"
	"os"

	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/gff"
	"github.com/biogo/store/interval"
)

func main() {
	flag.Usage = func() {
		fmt.Println(`usage: cull < infile.gff > outfile.gff`)
		os.Exit(0)
	}
	flag.Parse()
	r := gff.NewReader(os.Stdin)
	sc := featio.NewScanner(r)
	var feats []*gff.Feature
	for sc.Next() {
		feats = append(feats, sc.Feat().(*gff.Feature))
	}
	if err := sc.Error(); err != nil {
		log.Fatal(err)
	}
	w := gff.NewWriter(os.Stdout, 60, true)
	for _, f := range cullContained(feats) {
		_, err := w.Write(f)
		if err != nil {
			log.Fatal(err)
		}
	}
}

// cullContained returns a copy of hits with all hits that are completely contained by
// a higher scoring hit removed.
func cullContained(hits []*gff.Feature) []*gff.Feature {
	var tree interval.IntTree
	for i, f := range hits {
		if f.FeatScore == nil {
			continue
		}
		err := tree.Insert(subjectInterval{uid: uintptr(i), Feature: f}, true)
		if err != nil {
			log.Fatal(err)
		}
	}
	tree.AdjustRanges()
	var culled []*gff.Feature
outer:
	for _, f := range hits {
		if f.FeatScore != nil {
			o := tree.Get(subjectInterval{Feature: f})
			for _, h := range o {
				if *h.(subjectInterval).FeatScore > *f.FeatScore {
					continue outer
				}
			}
		}
		culled = append(culled, f)
	}
	return culled
}

type subjectInterval struct {
	uid uintptr
	*gff.Feature
}

// Overlap returns whether the b interval completely contains i.
func (i subjectInterval) Overlap(b interval.IntRange) bool {
	return b.Start <= i.FeatStart && i.FeatEnd <= b.End
}
func (i subjectInterval) ID() uintptr { return i.uid }
func (i subjectInterval) Range() interval.IntRange {
	return interval.IntRange{Start: i.FeatStart, End: i.FeatEnd}
}
