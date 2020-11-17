// Copyright ©2020 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Package blast provides types and functions for invoking NCBI+ BLAST
// and interpreting the returned results.
package blast

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"os/exec"
	"strconv"
	"strings"
	"text/template"

	"github.com/biogo/external"
)

type MakeDB struct {
	// Usage: makeblastdb -dbtype <type> -out <file>
	//
	// For details relating to options and parameters, see the BLAST manual.
	//
	Cmd string `buildarg:"{{if .}}{{.}}{{else}}makeblastdb{{end}}"` // blastn

	In          string `buildarg:"{{with .}}-in{{split}}{{.}}{{end}}"`                 // -in <s>
	Out         string `buildarg:"{{with .}}-out{{split}}{{.}}{{end}}"`                // -out <s>
	InputType   string `buildarg:"{{with .}}-input_type{{split}}{{.}}{{end}}"`         // -input_type <s>
	DBType      string `buildarg:"{{with .}}-dbtype{{split}}{{.}}{{end}}"`             // -dbtype <s>
	Title       string `buildarg:"{{with .}}-title{{split}}{{.}}{{end}}"`              // -title <s>
	ParseSeqids bool   `buildarg:"{{if .}}-parse_seqids{{end}}"`                       // -parse_seqids
	HashIndex   bool   `buildarg:"{{if .}}-hash_index{{end}}"`                         // -hash_index
	MaskData    string `buildarg:"{{with .}}-mask_data{{split}}{{.}}{{.}}{{end}}"`     // -mask_data <s>
	MaxFileSize string `buildarg:"{{with .}}-max_file_size{{split}}{{.}}{{.}}{{end}}"` // -max_file_size <s>
	TaxID       int    `buildarg:"{{with .}}-taxid{{split}}{{.}}{{end}}"`              // -taxid <n>
	TaxIDMap    string `buildarg:"{{with .}}-taxid_map{{split}}{{.}}{{end}}"`          // -taxid_map <s>
	LogFile     string `buildarg:"{{with .}}-logfile{{split}}{{.}}{{end}}"`            // -logfile <s>

	// ExtraFlags will be passed through to makeblastdb as flags.
	ExtraFlags string
}

func (m MakeDB) BuildCommand() (*exec.Cmd, error) {
	if m.DBType == "" {
		return nil, errors.New("makeblastdb: missing dbtype")
	}
	if m.Out == "" {
		return nil, errors.New("makeblastdb: missing out filename")
	}
	var extra []string
	if m.ExtraFlags != "" {
		extra = strings.Split(m.ExtraFlags, " ")
	}
	cl := external.Must(external.Build(m))
	return exec.Command(cl[0], append(cl[1:], extra...)...), nil
}

type Nucleic struct {
	// Usage: blastn -db <file> -query <file>
	//
	// For details relating to options and parameters, see the BLAST manual.
	//
	Cmd string `buildarg:"{{if .}}{{.}}{{else}}blastn{{end}}"` // blastn

	// Parameter:
	EValue        float64 `buildarg:"{{if .}}-evalue{{split}}{{.}}{{end}}"`          // -evalue <f.>
	WordSize      int     `buildarg:"{{if .}}-word_size{{split}}{{.}}{{end}}"`       // -word_size <n>
	Dust          *Dust   `buildarg:"{{if .}}-dust{{split}}{{dust .}}{{end}}"`       // -dust <...>
	SoftMask      bool    `buildarg:"-soft_masking{{split}}{{.}}"`                   // -soft_masking <b>
	Reward        int     `buildarg:"{{if .}}-reward{{split}}{{.}}{{end}}"`          // -reward <n>
	Penalty       int     `buildarg:"{{if .}}-penalty{{split}}{{.}}{{end}}"`         // -penalty <n>
	XdropUngap    int     `buildarg:"{{if .}}-xdrop_ungap{{split}}{{.}}{{end}}"`     // -xdrop_ungap <n>
	XdropGap      int     `buildarg:"{{if .}}-xdrop_gap{{split}}{{.}}{{end}}"`       // -xdrop_gap <n>
	XdropGapFinal int     `buildarg:"{{if .}}-xdrop_gap_final{{split}}{{.}}{{end}}"` // -xdrop_gap_final <n>
	GapOpen       int     `buildarg:"{{if .}}-gapopen{{split}}{{.}}{{end}}"`         // -gapopen <n>
	GapExtend     int     `buildarg:"{{if .}}-gapextend{{split}}{{.}}{{end}}"`       // -gapextend <n>
	NumAlignments int     `buildarg:"{{if .}}-num_alignments{{split}}{{.}}{{end}}"`  // -num_alignments <n>
	SearchSpace   int     `buildarg:"{{if .}}-searchsp{{split}}{{.}}{{end}}"`        // -searchsp <n>
	ParseDeflines bool    `buildarg:"{{if .}}-parse_deflines{{end}}"`                // -parse_deflines

	// Input:
	Query    string `buildarg:"-query{{split}}{{.}}"`                  // -query <s>
	Subject  string `buildarg:"{{if .}}-subject{{split}}{{.}}{{end}}"` // -subject <s>
	Database string `buildarg:"{{if .}}-db{{split}}{{.}}{{end}}"`      // -db <s>

	// Output:
	OutFormat int `buildarg:"{{if .}}-outfmt{{split}}{{.}}{{end}}"` // -outfmt <n>

	// Performance:
	Threads int `buildarg:"{{if .}}-num_threads{{split}}{{.}}{{end}}"` // -num_threads <n>

	// ExtraFlags will be passed through to blastn as flags.
	ExtraFlags string
}

func (n Nucleic) BuildCommand() (*exec.Cmd, error) {
	cl := external.Must(external.Build(n, template.FuncMap{"dust": dust}))
	var extra []string
	if n.ExtraFlags != "" {
		extra = strings.Split(n.ExtraFlags, " ")
	}
	return exec.Command(cl[0], append(cl[1:], extra...)...), nil
}

// Dust options.
type Dust struct {
	Filter bool
	Level  int
	Window int
	Linker int
}

func dust(d Dust) string {
	if !d.Filter {
		return "no"
	}
	if d.Level == 0 && d.Window == 0 && d.Linker == 0 {
		return "yes"
	}
	return fmt.Sprintf("%d %d %d", d.Level, d.Window, d.Linker)
}

type Record struct {
	QueryAccVer     string
	SubjectAccVer   string
	PctIdentity     float64
	AlignmentLength int
	Mismatches      int
	GapOpens        int
	QueryStart      int
	QueryEnd        int
	SubjectStart    int
	SubjectEnd      int
	EValue          float64
	BitScore        float64

	Strand int8

	// Iteration is the blast Iteration
	// that gave the blast hit.
	Iteration int `json:",omitempty"`

	// UID of hit connecting HSPs in a BLAST hit.
	UID int64 `json:",omitempty"`
}

func ParseTabular(r io.Reader, iteration int) ([]Record, error) {
	// column indices for default blast output tabular format 6 and 7.
	const (
		QueryAccVer = iota
		SubjectAccVer
		PctIdentity
		AlignmentLength
		Mismatches
		GapOpens
		QueryStart
		QueryEnd
		SubjectStart
		SubjectEnd
		EValue
		BitScore
		numFields
	)

	var recs []Record
	sc := bufio.NewScanner(r)
	for sc.Scan() {
		line := sc.Bytes()
		if bytes.HasPrefix(line, []byte("#")) {
			// Allow format 7 as well.
			continue
		}
		f := bytes.Split(line, []byte("\t"))
		if len(f) != numFields {
			return recs, fmt.Errorf("unexpected number of fields: %q", f)
		}

		// For some reason, NCBI think it's reasonable to sometimes
		// contaminate numeric fields with flanking whitespace.
		// So we trim whitespace from all fields just in case.
		r := Record{
			QueryAccVer:   string(bytes.TrimSpace(f[QueryAccVer])),
			SubjectAccVer: string(bytes.TrimSpace(f[SubjectAccVer])),
			Iteration:     iteration,
		}
		var err error
		r.PctIdentity, err = strconv.ParseFloat(string(bytes.TrimSpace(f[PctIdentity])), 64)
		if err != nil {
			return recs, fmt.Errorf("error in line: %s: %w", line, err)
		}
		r.AlignmentLength, err = strconv.Atoi(string(bytes.TrimSpace(f[AlignmentLength])))
		if err != nil {
			return recs, fmt.Errorf("error in line: %s: %w", line, err)
		}
		r.Mismatches, err = strconv.Atoi(string(bytes.TrimSpace(f[Mismatches])))
		if err != nil {
			return recs, fmt.Errorf("error in line: %s: %w", line, err)
		}
		r.GapOpens, err = strconv.Atoi(string(bytes.TrimSpace(f[GapOpens])))
		if err != nil {
			return recs, fmt.Errorf("error in line: %s: %w", line, err)
		}
		r.QueryStart, err = strconv.Atoi(string(bytes.TrimSpace(f[QueryStart])))
		if err != nil {
			return recs, fmt.Errorf("error in line: %s: %w", line, err)
		}
		r.QueryStart-- // Use zero-based indexing internally.
		r.QueryEnd, err = strconv.Atoi(string(bytes.TrimSpace(f[QueryEnd])))
		if err != nil {
			return recs, fmt.Errorf("error in line: %s: %w", line, err)
		}
		r.SubjectStart, err = strconv.Atoi(string(bytes.TrimSpace(f[SubjectStart])))
		if err != nil {
			return recs, fmt.Errorf("error in line: %s: %w", line, err)
		}
		r.SubjectStart-- // Use zero-based indexing internally.
		r.SubjectEnd, err = strconv.Atoi(string(bytes.TrimSpace(f[SubjectEnd])))
		if err != nil {
			return recs, fmt.Errorf("error in line: %s: %w", line, err)
		}
		r.EValue, err = strconv.ParseFloat(string(bytes.TrimSpace(f[EValue])), 64)
		if err != nil {
			return recs, fmt.Errorf("error in line: %s: %w", line, err)
		}
		r.BitScore, err = strconv.ParseFloat(string(bytes.TrimSpace(f[BitScore])), 64)
		if err != nil {
			return recs, fmt.Errorf("error in line: %s: %w", line, err)
		}
		r.Strand = 1
		if r.SubjectEnd < r.SubjectStart {
			r.Strand = -1
		}
		if r.QueryEnd < r.QueryStart {
			panic("inverted query")
		}
		recs = append(recs, r)
	}
	err := sc.Err()
	return recs, err
}
