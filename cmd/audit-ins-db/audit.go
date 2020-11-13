// Copyright ©2020 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// The audit-ins-db command allows the internal data stores generated during
// a run of ins to be queried. There are three persisted data stores in addition
// to the artifacts that are left by blastn.
//  - forward.db — the result of forward blast searches
//  - regions.db — aggregated hits obtained by merging closely located
//                 blast hits stored in forward.db
//  - reverse.db — the result of reciprocal blast searches
// The db files will be found in the working directory notes in the log output
// of ins and will remain after ins completes an analysis if it is given the
// -work flag.
// Each of the databases must be named as described here for audit-ins-db to
// understand their contents. Output from audit-ins-db is a JSON stream on stdout.
//
// forward.db and reverse.db
//
// The forward.db and reverse.db files contains BLAST hit results in JSON
// corresponding to the following Go struct. The query fields refer to the
// identified repeat family and the subject fields refer to the identified
// genomic region.
//  struct {
//  	SubjectAccVer string
//  	SubjectLeft   int64
//  	SubjectRight  int64
//  	QueryAccVer   string
//  	QueryStart    int64
//  	QueryEnd      int64
//  	BitScore      float64
//  	Strand        int8
//  }
//
// regions.db
//
// The regions.db file contains aggregated blast hit results in JSON
// corresponding to the following Go struct. The query field refers to the
// identified repeat family and the subject fields refer to the identified
// genomic region. Count indicates the number of hits that were included in
// the region.
//  struct {
//  	SubjectAccVer string
//  	SubjectLeft   int64
//  	SubjectRight  int64
//  	QueryAccVer   string
//  	Strand        int8
//  	Count         int64
//  }
package main

import (
	"encoding/binary"
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"

	"modernc.org/kv"

	"github.com/kortschak/ins/internal/store"
)

var order = binary.BigEndian

func main() {
	path := flag.String("db", "", "specify db file to audit (base must match '{forward,regions,reverse}.db')")
	flag.Parse()
	base := filepath.Base(*path)
	switch base {
	case "forward.db", "regions.db", "reverse.db":
	default:
		flag.Usage()
		os.Exit(2)
	}

	var enc *json.Encoder
	if base == "regions.db" {
		enc = json.NewEncoder(os.Stdout)
	}

	orderFor := map[string]func(x, y []byte) int{
		"forward.db": store.GroupByQueryOrderSubjectLeft,
		"regions.db": store.GroupByQueryOrderSubjectLeft,
		"reverse.db": store.BySubjectPosition,
	}
	opts := &kv.Options{Compare: orderFor[base]}
	db, err := kv.Open(*path, opts)
	if err != nil {
		log.Fatal(err)
	}
	defer db.Close()

	it, err := db.SeekFirst()
	if err != nil {
		if err == io.EOF {
			return
		}
		log.Fatal(err)
	}
	for {
		k, v, err := it.Next()
		if err != nil {
			if err == io.EOF {
				break
			}
			log.Fatal(err)
		}
		switch base {
		case "forward.db", "reverse.db":
			os.Stdout.Write(v)
			fmt.Println()
		case "regions.db":
			r := store.UnmarshalBlastRecordKey(k)
			n := int64(order.Uint64(v))
			err = enc.Encode(region{
				SubjectAccVer: r.SubjectAccVer,
				SubjectLeft:   r.SubjectLeft,
				SubjectRight:  r.SubjectRight,
				QueryAccVer:   r.QueryAccVer,
				Strand:        r.Strand,
				Count:         n,
			})
			if err != nil {
				log.Fatal(err)
			}
		default:
			panic("unreachable")
		}
	}
}

type region struct {
	SubjectAccVer string
	SubjectLeft   int64
	SubjectRight  int64
	QueryAccVer   string
	Strand        int8
	Count         int64
}
