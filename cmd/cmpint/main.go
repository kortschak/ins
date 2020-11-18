// Copyright ©2020 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// The cmpint program compares the intervals in two files. It takes two
// GTF file inputs describing repeat annotations and compares overlapping
// regions. The output of the analysis is the number of bases that agree
// between the inputs, the number of bases that are covered in one, but
// not the other, and the number of bases where the annotation differs.
// These analyses are done for both the repeat type and the repeat class,
// and are emitted on stdout as a JSON object.
//
// If a dot flag is provided, descriptions of the discordances between the
// feature sets as a graph in DOT format, with edge weights representing
// counts of mismatched bases.
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"sort"
	"strings"

	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/gff"
	"github.com/biogo/store/step"
	"gonum.org/v1/gonum/graph"
	"gonum.org/v1/gonum/graph/encoding"
	"gonum.org/v1/gonum/graph/encoding/dot"
	"gonum.org/v1/gonum/graph/simple"
)

func main() {
	aFile := flag.String("a", "", "specify the input file a name (required)")
	bFile := flag.String("b", "", "specify the input file b name (required)")
	out := flag.String("dot", "", "specify prefix for DOT files describing disagreements")
	none := flag.String("none", "none", "specify label for 'no annotation")

	flag.Parse()
	if *aFile == "" || *bFile == "" {
		flag.Usage()
		os.Exit(2)
	}

	chrs := make(map[string]bool)
	types := make(map[string]*step.Vector)
	classes := make(map[string]*step.Vector)
	err := steps(*aFile, func(f *gff.Feature) error {
		chrs[f.SeqName] = true

		typ, class, err := typeClassOf(f)
		if err != nil {
			return err
		}

		tv, ok := types[f.SeqName]
		if !ok {
			tv, err = step.New(0, 1, pair{})
			if err != nil {
				return err
			}
			tv.Relaxed = true
			types[f.SeqName] = tv
		}
		err = tv.ApplyRange(f.FeatStart, f.FeatEnd, func(e step.Equaler) step.Equaler {
			t := e.(pair)
			if f.FeatScore != nil && *f.FeatScore > t.aScore {
				t.a = typ
				t.aScore = *f.FeatScore
			}
			return t
		})
		if err != nil {
			return err
		}

		cv, ok := classes[f.SeqName]
		if !ok {
			cv, err = step.New(0, 1, pair{})
			if err != nil {
				return err
			}
			cv.Relaxed = true
			classes[f.SeqName] = cv
		}
		err = cv.ApplyRange(f.FeatStart, f.FeatEnd, func(e step.Equaler) step.Equaler {
			c := e.(pair)
			if f.FeatScore != nil && *f.FeatScore > c.aScore {
				c.a = class
				c.aScore = *f.FeatScore
			}
			return c
		})
		return err
	})
	if err != nil {
		log.Fatal(err)
	}
	err = steps(*bFile, func(f *gff.Feature) error {
		chrs[f.SeqName] = true

		typ, class, err := typeClassOf(f)
		if err != nil {
			return err
		}

		tv, ok := types[f.SeqName]
		if !ok {
			tv, err = step.New(0, 1, pair{})
			if err != nil {
				return err
			}
			tv.Relaxed = true
			types[f.SeqName] = tv
		}
		err = tv.ApplyRange(f.FeatStart, f.FeatEnd, func(e step.Equaler) step.Equaler {
			t := e.(pair)
			if f.FeatScore != nil && *f.FeatScore > t.bScore {
				t.b = typ
				t.bScore = *f.FeatScore
			}
			return t
		})
		if err != nil {
			return err
		}

		cv, ok := classes[f.SeqName]
		if !ok {
			cv, err = step.New(0, 1, pair{})
			if err != nil {
				return err
			}
			cv.Relaxed = true
			classes[f.SeqName] = cv
		}
		err = cv.ApplyRange(f.FeatStart, f.FeatEnd, func(e step.Equaler) step.Equaler {
			c := e.(pair)
			if f.FeatScore != nil && *f.FeatScore > c.bScore {
				c.b = class
				c.bScore = *f.FeatScore
			}
			return c
		})
		return err
	})
	if err != nil {
		log.Fatal(err)
	}

	// FIXME: do we need this; currently not
	var chroms []string
	for c := range chrs {
		chroms = append(chroms, c)
	}
	sort.Strings(chroms)

	// Process classes.
	var (
		classAgree      int
		aMissingClass   int
		bMissingClass   int
		classMismatch   int
		classMismatches = make(map[names]int)
	)
	for _, chr := range chroms {
		classes[chr].Do(func(start, end int, e step.Equaler) {
			c := e.(pair)
			if c.isZero() {
				return
			}
			len := end - start
			switch {
			case c.a == c.b:
				classAgree += len
			case c.a == "":
				aMissingClass += len
				classMismatches[names{a: "", b: c.b}] += len
			case c.b == "":
				bMissingClass += len
				classMismatches[names{a: c.a, b: ""}] += len
			default:
				classMismatch += len
				classMismatches[c.names] += len
			}
		})
	}

	// Process types.
	var (
		typeAgree      int
		aMissingType   int
		bMissingType   int
		typeMismatch   int
		typeMismatches = make(map[names]int)
	)
	for _, chr := range chroms {
		types[chr].Do(func(start, end int, e step.Equaler) {
			t := e.(pair)
			if t.isZero() {
				return
			}
			len := end - start
			switch {
			case t.a == t.b:
				typeAgree += len
			case t.a == "":
				aMissingType += len
				typeMismatches[names{a: "", b: t.b}] += len
			case t.b == "":
				bMissingType += len
				typeMismatches[names{a: t.a, b: ""}] += len
			default:
				typeMismatch += len
				typeMismatches[t.names] += len
			}
		})
	}

	type record struct {
		Agree    int `json:"agree"`
		AMissing int `json:"a-missing"`
		BMissing int `json:"b-missing"`
		Mismatch int `json:"mismatch"`
	}
	type report struct {
		Class record `json:"class"`
		Type  record `json:"type"`
	}

	m, err := json.Marshal(report{
		Class: record{
			Agree:    classAgree,
			AMissing: aMissingClass,
			BMissing: bMissingClass,
			Mismatch: classMismatch,
		},
		Type: record{
			Agree:    typeAgree,
			AMissing: aMissingType,
			BMissing: bMissingType,
			Mismatch: typeMismatch,
		},
	})
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("%s\n", m)
	if *out != "" {
		err = dotOut(*out+".class.dot", *aFile, *bFile, classMismatches, *none)
		if err != nil {
			log.Fatal(err)
		}
		err = dotOut(*out+".type.dot", *aFile, *bFile, typeMismatches, *none)
		if err != nil {
			log.Fatal(err)
		}
	}
}

func steps(path string, fn func(*gff.Feature) error) error {
	f, err := os.Open(path)
	if err != nil {
		return err
	}
	r := gff.NewReader(f)
	sc := featio.NewScanner(r)
	for sc.Next() {
		gf := sc.Feat().(*gff.Feature)
		err := fn(gf)
		if err != nil {
			return err
		}
	}
	return sc.Error()
}

func typeClassOf(f *gff.Feature) (typ, class string, err error) {
	attr := f.FeatAttributes.Get("Repeat")
	fields := strings.Split(attr, " ")
	if len(fields) != 5 {
		return "", "", fmt.Errorf("unexpected repeat attribute syntax: %q", attr)
	}
	return fields[0], fields[1], nil
}

// pair is a step vector element with two string values
// which are either repeat element type or class.
type pair struct {
	names

	aScore float64
	bScore float64
}

type names struct {
	a, b string
}

func (p pair) isZero() bool {
	return p.names == names{}
}

func (p pair) Equal(e step.Equaler) bool {
	return p.names == e.(pair).names
}

func dotOut(path, aFile, bFile string, edges map[names]int, none string) error {
	g := newNameGraph(none)
	for p, w := range edges {
		e := edge{
			f: g.nodeFor(aFile, p.a),
			t: g.nodeFor(bFile, p.b),
			w: float64(w),
		}
		g.SetWeightedEdge(e)
	}
	b, err := dot.Marshal(g, "discord", "", "\t")
	if err != nil {
		return err
	}
	err = ioutil.WriteFile(path, b, 0o664)
	if err != nil {
		return err
	}
	return nil
}

type nameGraph struct {
	*simple.WeightedUndirectedGraph
	idFor map[string]int64
	none  string
}

func newNameGraph(none string) nameGraph {
	return nameGraph{
		WeightedUndirectedGraph: simple.NewWeightedUndirectedGraph(0, 0),
		idFor:                   make(map[string]int64),
		none:                    none,
	}
}

func (g nameGraph) nodeFor(file, s string) graph.Node {
	if s == "" {
		s = g.none
	}
	s = file + ":" + s
	id, ok := g.idFor[s]
	if ok {
		return g.Node(id)
	}
	id = g.WeightedUndirectedGraph.NewNode().ID()
	g.idFor[s] = id
	n := node{id: id, name: s}
	g.AddNode(n)
	return n
}

type node struct {
	id   int64
	name string
}

func (n node) ID() int64     { return n.id }
func (n node) DOTID() string { return n.name }

type edge struct {
	f, t graph.Node
	w    float64
}

func (e edge) From() graph.Node         { return e.f }
func (e edge) To() graph.Node           { return e.t }
func (e edge) ReversedEdge() graph.Edge { return edge{f: e.t, t: e.f, w: e.w} }
func (e edge) Weight() float64          { return e.w }
func (e edge) Attributes() []encoding.Attribute {
	return []encoding.Attribute{{Key: "weight", Value: fmt.Sprint(e.w)}}
}
