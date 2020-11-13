// Copyright ©2020 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package store

import (
	"bytes"
	"encoding/binary"
	"math"

	"github.com/kortschak/ins/blast"
)

// GroupByQueryOrderSubjectLeft is a kv compare function, ordering by strand, query name,
// subject name, subject position and BLAST bitscore.
func GroupByQueryOrderSubjectLeft(x, y []byte) int {
	if bytes.Equal(x, y) {
		return 0
	}

	rx := UnmarshalBlastRecordKey(x)
	ry := UnmarshalBlastRecordKey(y)

	// Separate strands, (+) first.
	switch {
	case rx.Strand > ry.Strand:
		return -1
	case rx.Strand < ry.Strand:
		return 1
	}

	// Group elements of the same type.
	switch {
	case rx.QueryAccVer < ry.QueryAccVer:
		return -1
	case rx.QueryAccVer > ry.QueryAccVer:
		return 1
	}

	// Sort by left position, with higher scoring matches first.
	switch {
	case rx.SubjectAccVer < ry.SubjectAccVer:
		return -1
	case rx.SubjectAccVer > ry.SubjectAccVer:
		return 1
	}
	switch {
	case rx.SubjectLeft < ry.SubjectLeft:
		return -1
	case rx.SubjectLeft > ry.SubjectLeft:
		return 1
	}
	switch {
	case rx.SubjectRight < ry.SubjectRight:
		return -1
	case rx.SubjectRight > ry.SubjectRight:
		return 1
	}
	switch {
	case rx.BitScore > ry.BitScore:
		return -1
	case rx.BitScore < ry.BitScore:
		return 1
	}

	// Ensure key uniqueness.
	switch {
	case rx.QueryStart < ry.QueryStart:
		return -1
	case rx.QueryStart > ry.QueryStart:
		return 1
	}
	switch {
	case rx.QueryEnd < ry.QueryEnd:
		return -1
	case rx.QueryEnd > ry.QueryEnd:
		return 1
	}

	panic("unreachable")
}

// BySubjectPosition is a kv compare function, ordering by strand, subject name,
// subject position and BLAST bitscore.
func BySubjectPosition(x, y []byte) int {
	if bytes.Equal(x, y) {
		return 0
	}

	rx := UnmarshalBlastRecordKey(x)
	ry := UnmarshalBlastRecordKey(y)

	// Separate strands, (+) first.
	switch {
	case rx.Strand > ry.Strand:
		return -1
	case rx.Strand < ry.Strand:
		return 1
	}

	// Sort by left position, longer repeats first,
	// and with higher scoring matches first.
	switch {
	case rx.SubjectAccVer < ry.SubjectAccVer:
		return -1
	case rx.SubjectAccVer > ry.SubjectAccVer:
		return 1
	}
	switch {
	case rx.SubjectLeft < ry.SubjectLeft:
		return -1
	case rx.SubjectLeft > ry.SubjectLeft:
		return 1
	}
	switch {
	case rx.SubjectRight > ry.SubjectRight:
		return -1
	case rx.SubjectRight < ry.SubjectRight:
		return 1
	}
	switch {
	case rx.BitScore > ry.BitScore:
		return -1
	case rx.BitScore < ry.BitScore:
		return 1
	}

	// Ensure key uniqueness.
	switch {
	case rx.QueryAccVer < ry.QueryAccVer:
		return -1
	case rx.QueryAccVer > ry.QueryAccVer:
		return 1
	}
	switch {
	case rx.QueryStart < ry.QueryStart:
		return -1
	case rx.QueryStart > ry.QueryStart:
		return 1
	}
	switch {
	case rx.QueryEnd < ry.QueryEnd:
		return -1
	case rx.QueryEnd > ry.QueryEnd:
		return 1
	}

	panic("unreachable")
}

// MarshalInt returns a slice encoding n as an int64.
func MarshalInt(n int) []byte {
	var buf [8]byte
	order.PutUint64(buf[:], uint64(n))
	return buf[:]
}

type BlastRecordKey struct {
	SubjectAccVer string
	SubjectLeft   int64
	SubjectRight  int64
	QueryAccVer   string
	QueryStart    int64
	QueryEnd      int64
	BitScore      float64
	Strand        int8
}

var order = binary.BigEndian

func MarshalBlastRecordKey(r blast.Record) []byte {
	var (
		buf bytes.Buffer
		b   [8]byte
	)
	order.PutUint64(b[:], uint64(len(r.SubjectAccVer)))
	buf.Write(b[:])
	buf.WriteString(r.SubjectAccVer)
	left := r.SubjectStart
	right := r.SubjectEnd
	if right < left {
		left, right = right, left
		// The query coordinates are used for disambiguation.
		// This reordering retains the relative coordinates
		// for cases where there are matches in both orientations.
		r.QueryStart, r.QueryEnd = r.QueryEnd, r.QueryStart
	}
	order.PutUint64(b[:], uint64(left))
	buf.Write(b[:])
	order.PutUint64(b[:], uint64(right))
	buf.Write(b[:])

	order.PutUint64(b[:], uint64(len(r.QueryAccVer)))
	buf.Write(b[:])
	buf.WriteString(r.QueryAccVer)
	order.PutUint64(b[:], uint64(r.QueryStart))
	buf.Write(b[:])
	order.PutUint64(b[:], uint64(r.QueryEnd))
	buf.Write(b[:])
	order.PutUint64(b[:], math.Float64bits(r.BitScore))
	buf.Write(b[:])
	buf.WriteByte(byte(r.Strand))
	return buf.Bytes()
}

func UnmarshalBlastRecordKey(data []byte) BlastRecordKey {
	var k BlastRecordKey
	n64 := binary.Size(uint64(0))
	n := order.Uint64(data[:n64])
	data = data[n64:]
	k.SubjectAccVer = string(data[:n])
	data = data[n:]
	k.SubjectLeft = int64(order.Uint64(data[:n64]))
	data = data[n64:]
	k.SubjectRight = int64(order.Uint64(data[:n64]))
	data = data[n64:]
	n = order.Uint64(data[:n64])
	data = data[n64:]
	k.QueryAccVer = string(data[:n])
	data = data[n:]
	k.QueryStart = int64(order.Uint64(data[:n64]))
	data = data[n64:]
	k.QueryEnd = int64(order.Uint64(data[:n64]))
	data = data[n64:]
	k.BitScore = math.Float64frombits(order.Uint64(data[:n64]))
	data = data[n64:]
	k.Strand = int8(data[0])
	return k
}
