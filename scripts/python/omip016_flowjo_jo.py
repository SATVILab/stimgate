#!/usr/bin/env python3
"""Extract manual gates and compensation from the OMIP-016 FlowJo workspace.

The OMIP-016 workspace (FlowRepository FR-FCM-ZZ2T,
``PHY1002_20110118_ICS_SG_OMIP.jo``) is a legacy binary FlowJo for Mac
(v8-era, big-endian) workspace. CytoML and FlowKit read only XML workspaces,
so this script decodes the parts needed to re-apply the manual gating:

* the named compensation (spillover) matrix the samples use;
* for each requested sample, the population tree, each gate's parameters and
  polygon vertices (float32, big-endian, in the gate's data scale), and each
  Boolean gate's expression and references.

The decoding rules were derived by inspecting this one file. Every rule is
checked as it is applied (tree child counts must be consumed exactly, every
population must have exactly one gate definition, Boolean references must
resolve to siblings, the matrix must be square with a unit diagonal) and the
script fails rather than guess. It is not a general FlowJo reader.

Only the Python standard library is used. Usage::

    python3 -I omip016_flowjo_jo.py WORKSPACE.jo OUT_DIR SAMPLE.fcs [SAMPLE.fcs ...]

Outputs in OUT_DIR: ``compensation.csv``, ``populations.csv``,
``vertices.csv`` and ``workspace.json``.
"""

import csv
import hashlib
import json
import os
import re
import struct
import sys

# Pre-order population record: the 4-byte count of the parent's children
# (non-zero only on the first child), two 4-byte fields, then the Pascal
# name and the Pascal owner string.
POP_RE = re.compile(
    rb"\x00\x00\x00([\x00-\xff])\x00\x00\x00\x03\x00\x00\x00[\x00-\xff]"
    rb"([\x01-\x3f])([ -~]+?)([\x01-\x3f])([ -~]+?)\x00\x00\x00\x00\x00\x00\x00"
)
GEOM_TAGS = (b"Tpol", b"Tpy4")
FLOAT_LIMIT = 1e7


def u32(b, o):
    return struct.unpack(">I", b[o:o + 4])[0]


def read_lp_string(b, o):
    n = u32(b, o)
    return b[o + 4:o + 4 + n].decode("latin-1"), o + 4 + n


def compensation(b):
    """Return (name, channels, matrix) of the text spillover matrix."""
    m = re.search(rb"(mtx [ -~]+?)\r<\t>\r", b)
    if m is None:
        raise ValueError("Compensation matrix definition not found.")
    name = m.group(1).decode("latin-1")
    o = m.end()
    end = b.index(b"\x00", o)
    lines = b[o:end].decode("latin-1").split("\r")
    channels = lines[0].split("\t")
    n = len(channels)
    rows = [[float(v) for v in line.split("\t")] for line in lines[1:n + 1]]
    if len(rows) != n or any(len(r) != n for r in rows):
        raise ValueError("Compensation matrix is not square.")
    if any(abs(rows[i][i] - 1) > 1e-12 for i in range(n)):
        raise ValueError("Compensation matrix diagonal is not one.")
    return name, channels, rows


def sample_blocks(b, samples):
    """Byte range of each sample's population tree."""
    starts = sorted(
        (m.start(), m.group(1).decode("latin-1"))
        for m in re.finditer(rb"ICS ON FRESH CELLS/([ -~]+?\.fcs)", b)
    )
    out = {}
    for s in samples:
        idx = [i for i, (_, f) in enumerate(starts) if f == s]
        if len(idx) != 1:
            raise ValueError("Expected one workspace entry for %s; found %d." % (s, len(idx)))
        i = idx[0]
        lo = starts[i][0]
        # The sample's own trailing keyword record (its file name again)
        # closes the tree.
        hi = b.find(s.encode("latin-1"), lo + len(s) + 20)
        if hi < 0:
            hi = starts[i + 1][0] if i + 1 < len(starts) else len(b)
        out[s] = (lo, hi)
    return out


def parse_geometry(b, lo, hi, name):
    for tag in GEOM_TAGS:
        o = b.find(tag + name.encode("latin-1") + b"\r", lo, hi)
        if o >= 0:
            break
    else:
        return None
    p = b.find(b"\x01\x01\x00\x00", o, o + 80)
    if p < 0:
        raise ValueError("Gate '%s': parameter header not found." % name)
    q = b.index(b"\r", p + 4)
    params = b[p + 4:q].decode("latin-1").split("\t")
    if len(params) != 2:
        raise ValueError("Gate '%s': expected two parameters, got %r." % (name, params))
    n = u32(b, q + 1)
    if not 3 <= n <= 100:
        raise ValueError("Gate '%s': implausible vertex count %d." % (name, n))
    vals = struct.unpack(">%df" % (2 * n), b[q + 5:q + 5 + 8 * n])
    if any(not abs(v) < FLOAT_LIMIT for v in vals):
        raise ValueError("Gate '%s': implausible vertex values." % name)
    return {
        "type": "polygon",
        "tag": tag.decode(),
        "params": params,
        "vertices": [(vals[2 * i], vals[2 * i + 1]) for i in range(n)],
        "offset": o,
    }


def parse_boolean(b, lo, hi):
    o = b.find(b"Bool>", lo, hi)
    if o < 0:
        return None
    # Expression length follows the header floats; find it by matching a
    # length-prefixed string that starts with "G".
    m = re.compile(rb"\x00\x00\x00([\x01-\xff])(G[0-9 |&!()G]+)").search(b, o, o + 80)
    if m is None or len(m.group(2)) != m.group(1)[0]:
        raise ValueError("Boolean expression not decoded at %d." % o)
    expr = m.group(2).decode("latin-1")
    p = m.end()
    nref = u32(b, p)
    p += 4
    refs = []
    for _ in range(nref):
        s, p = read_lp_string(b, p)
        refs.append(s)
    used = sorted({int(g) for g in re.findall(r"G(\d+)", expr)})
    if used and max(used) >= nref:
        raise ValueError("Boolean '%s' refers to a missing operand." % expr)
    return {"type": "boolean", "expr": expr, "refs": refs, "offset": o}


def parse_sample(b, lo, hi):
    recs = []
    for m in POP_RE.finditer(b, lo, hi):
        if m.group(2)[0] != len(m.group(3)) or m.group(4)[0] != len(m.group(5)):
            continue
        recs.append({
            "count": m.group(1)[0],
            "name": m.group(3).decode("latin-1"),
            "start": m.start(),
        })
    if not recs:
        raise ValueError("No populations found.")
    # Rebuild the tree from the pre-order child counts: a non-zero count on
    # a record means it is the first of that many children of the record
    # before it (of the sample root, id 0, for the first record).
    if recs[0]["count"] < 1:
        raise ValueError("First population lacks a child count.")
    stack = []
    pops = []
    for i, r in enumerate(recs):
        if r["count"] > 0:
            stack.append({"id": i, "left": r["count"]})
        while stack and stack[-1]["left"] == 0:
            stack.pop()
        if not stack:
            raise ValueError("Population '%s' outside the tree." % r["name"])
        stack[-1]["left"] -= 1
        pops.append({"id": i + 1, "parent": stack[-1]["id"], "name": r["name"]})
    if any(s["left"] for s in stack):
        raise ValueError("Population tree child counts were not consumed.")

    by_id = {p["id"]: p for p in pops}

    def path(p):
        parts = []
        while p["id"] != 0:
            parts.append(p["name"])
            p = by_id.get(p["parent"], {"id": 0})
        return "/".join(reversed(parts))

    for i, (p, r) in enumerate(zip(pops, recs)):
        rlo = r["start"]
        rhi = recs[i + 1]["start"] if i + 1 < len(recs) else hi
        geom = parse_geometry(b, rlo, rhi, p["name"])
        boolean = parse_boolean(b, rlo, rhi)
        if (geom is None) == (boolean is None):
            raise ValueError("Population '%s' needs exactly one gate definition." % p["name"])
        p.update(geom or boolean)
        p["path"] = path(p)
    for p in pops:
        if p["type"] == "boolean":
            siblings = {q["name"] for q in pops if q["parent"] == p["parent"]}
            for ref in p["refs"]:
                if not ref.startswith("/") or ref[1:] not in siblings:
                    raise ValueError("Boolean reference '%s' of '%s' is not a sibling." % (ref, p["name"]))
    return pops


def main(argv):
    if len(argv) < 4:
        sys.exit(__doc__)
    path_jo, out_dir, samples = argv[1], argv[2], argv[3:]
    with open(path_jo, "rb") as f:
        b = f.read()
    if not b.startswith(b"FlowJo"):
        raise ValueError("Not a binary FlowJo workspace: %s" % path_jo)
    comp_name, comp_chnl, comp = compensation(b)
    blocks = sample_blocks(b, samples)
    os.makedirs(out_dir, exist_ok=True)

    with open(os.path.join(out_dir, "compensation.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["channel"] + comp_chnl)
        for ch, row in zip(comp_chnl, comp):
            w.writerow([ch] + [repr(v) for v in row])

    pop_rows, vert_rows, sample_meta = [], [], []
    for s in samples:
        lo, hi = blocks[s]
        # The sample's keyword records (just before and inside its block)
        # name the compensation matrix it uses; each must be the decoded one.
        refs = [r.decode("latin-1") for r in re.findall(
            rb"FJ_CompMatrixName[^ -~]?([ -~]{%d})" % len(comp_name),
            b[max(0, lo - 4000):hi],
        )]
        if not refs or any(r != comp_name for r in refs):
            raise ValueError("Sample %s does not use compensation '%s': %r" % (s, comp_name, refs))
        sample_meta.append({
            "sample": s,
            "byteRange": [lo, hi],
            "compMatrixReferences": len(refs),
        })
        for p in parse_sample(b, lo, hi):
            pop_rows.append({
                "sample": s,
                "pop_id": p["id"],
                "parent_id": p["parent"],
                "name": p["name"],
                "path": p["path"],
                "type": p["type"],
                "x_param": p["params"][0] if p["type"] == "polygon" else "",
                "y_param": p["params"][1] if p["type"] == "polygon" else "",
                "n_vertices": len(p["vertices"]) if p["type"] == "polygon" else 0,
                "bool_expr": p.get("expr", ""),
                "bool_refs": ";".join(p.get("refs", [])),
                "byte_offset": p["offset"],
            })
            for k, (x, y) in enumerate(p.get("vertices", []), start=1):
                vert_rows.append({
                    "sample": s, "pop_id": p["id"], "path": p["path"],
                    "vertex": k, "x": repr(x), "y": repr(y),
                })

    for name, rows in (("populations.csv", pop_rows), ("vertices.csv", vert_rows)):
        with open(os.path.join(out_dir, name), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

    meta = {
        "workspace": os.path.basename(path_jo),
        "workspaceMd5": hashlib.md5(b).hexdigest(),
        "format": "FlowJo for Mac binary workspace (big-endian)",
        "parser": os.path.basename(__file__),
        "compensation": {"name": comp_name, "channels": comp_chnl},
        "samples": sample_meta,
    }
    with open(os.path.join(out_dir, "workspace.json"), "w") as f:
        json.dump(meta, f, indent=2)


if __name__ == "__main__":
    main(sys.argv)
