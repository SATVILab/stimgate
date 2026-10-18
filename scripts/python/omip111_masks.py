"""Reconstruct OMIP-111 parent and cytokine gate membership using FlowKit 1.3.2.

Export original FCS event indices, never the clipped FlowJo display values.
Audit every stored gate count. The R analysis applies a common raw-data transform.

The authors' cytokine gates are cytokine-by-CD44 rectangles inside a non-naive
(CD44 versus CD62L) gate. The CD44+ parent of each cytokine follows their
gating path exactly: the non-naive gate restricted to that tube's cytokine
gate's CD44 interval, compared on the gate's display scale. Exported rows are
the non-naive cells, with one parent column per cytokine. The same CD44
interval within the whole lineage gate (ignoring CD62L) is counted for
comparison, because stimulation sheds CD62L.
"""
import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import flowkit as fk
import numpy as np
import pandas as pd

CD44 = "[RealBlue 744]-A"
CHANNELS = {"IFNg": "AF488-A", "IL2": "PE-A", "TNF": "PE-Cy7-A", "IL4_5": "BV421-A", "IL17A": "R718-A"}
GATES = {"IFNg": "IFNg+", "IL2": "IL-2", "TNF": "TNF+",
         "IL4_5": "Th2 Cells", "IL17A": "Th17 Cells"}
MARKERS = {"CD4": list(GATES), "CD8": list(GATES)[:3]}


def stored_gates(node, path=("root",), definitions=None):
    result = {}
    for child in node.findall("./Subpopulations/Population"):
        name = child.attrib["name"]
        result[(name, path)] = int(child.attrib["count"])
        if definitions is not None:
            definitions[(name, path)] = child
        result.update(stored_gates(child, path + (name,), definitions))
    return result


def raw_limits(definition):
    """Raw-unit rectangle limits stored in the FlowJo XML, by channel."""
    limits = {}
    for dimension in [d for d in definition.iter() if d.tag.rsplit("}", 1)[-1] == "dimension"]:
        channel_nodes = [n for n in dimension.iter() if n.tag.rsplit("}", 1)[-1] == "fcs-dimension"]
        if len(channel_nodes) != 1:
            raise ValueError("Unexpected author gate dimension")
        channel = next(v for k, v in channel_nodes[0].attrib.items() if k.rsplit("}", 1)[-1] == "name")
        limits[channel] = {k.rsplit("}", 1)[-1]: float(v) for k, v in dimension.attrib.items()}
    return limits


def export(raw_dir, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    audits, parents, references = [], [], []
    for strain, workspace in {"C57": "OMIP-ICS-C57BL_6.wsp", "BALB": "OMIP-ICS-BALB_c.wsp"}.items():
        xml = ET.parse(raw_dir / workspace)
        nodes = {n.attrib["name"]: n for n in xml.findall(".//Sample/SampleNode")}
        files = sorted(raw_dir.glob(f"[EF]* {strain}_M*_*.fcs"))
        wsp = fk.Workspace(str(raw_dir / workspace), fcs_samples=[str(f) for f in files])
        if len(files) != 10 or set(wsp.get_sample_ids()) != set(nodes):
            raise ValueError("Workspace/full-panel sample matching failed")
        for sid in sorted(nodes):
            if wsp.get_sample(sid).event_count != int(nodes[sid].attrib["count"]):
                raise ValueError("FCS/workspace event count mismatch")
            if wsp.get_comp_matrix(sid) is not None:
                raise ValueError("Unexpected compensation on already-unmixed FCS")
            wsp.analyze_samples(sample_id=sid, use_mp=False)
            definitions = {}
            stored = stored_gates(nodes[sid], definitions=definitions)
            sample = wsp.get_sample(sid)
            raw_events = sample.get_events(source="raw")
            transforms = wsp.get_transforms(sid)
            for name, path in wsp.get_gate_ids(sid):
                observed = int(wsp.get_gate_membership(sid, name, gate_path=path).sum())
                expected = stored[(name, path)]
                audits.append(dict(sample=sid, pop="/" + "/".join((path + (name,))[1:]),
                                   count_import=observed, count_flowjo=expected, difference=observed-expected))
            cd44_raw = raw_events[:, sample.get_channel_index(CD44)]
            for population in ("CD4", "CD8"):
                paths = [p for p in wsp.find_matching_gate_paths(sid, "Non-naive")
                         if p[-1] == f"{population}+ T Cells"]
                if len(paths) != 1:
                    raise ValueError("Ambiguous non-naive parent")
                path = paths[0]
                child_path = path + ("Non-naive",)
                lineage = wsp.get_gate_membership(sid, path[-1], gate_path=path[:-1])
                non_naive = wsp.get_gate_membership(sid, "Non-naive", gate_path=path)
                markers = MARKERS[population]
                gates = {m: wsp.get_gate(sid, GATES[m], gate_path=child_path) for m in markers}
                cd44_dims = {}
                for marker, gate in gates.items():
                    dims = {d.id: d for d in gate.dimensions}
                    if set(dims) != {CHANNELS[marker], CD44}:
                        raise ValueError("Author cytokine gate is not cytokine-by-CD44")
                    cd44_dims[marker] = dims[CD44]
                mask = non_naive
                table = pd.DataFrame({"eventIndex": np.flatnonzero(mask) + 1,
                                      "nonNaive": non_naive[mask].astype(int)})
                filename = sid.replace(".fcs", "") + f"-{population}.csv"
                for marker in markers:
                    full = wsp.get_gate_membership(sid, GATES[marker], gate_path=child_path)
                    dim = cd44_dims[marker]
                    display = transforms[dim.transformation_ref].apply(cd44_raw)
                    in_interval = (display >= dim.min) & (display < dim.max)
                    in_parent = non_naive & in_interval
                    if np.any(full & ~in_parent):
                        raise ValueError("Author cytokine-positive cells outside the CD44 parent")
                    limits = raw_limits(definitions[(GATES[marker], child_path)])
                    cytokine = raw_events[mask, sample.get_channel_index(CHANNELS[marker])]
                    cd44 = cd44_raw[mask]
                    parent = in_parent[mask]
                    lower, upper = limits[CHANNELS[marker]]["min"], limits[CHANNELS[marker]]["max"]
                    projected = parent & (cytokine > lower)
                    # Raw-unit CD44 interval and rectangle, for an independent
                    # check in R; they differ from the display-scale gates only
                    # at the edges.
                    raw_parent = non_naive[mask] & (cd44 >= limits[CD44]["min"]) & (cd44 < limits[CD44]["max"])
                    raw_rectangle = (non_naive[mask] & (cytokine >= lower) & (cytokine < upper) & raw_parent)
                    author = full[mask]
                    table[f"parent_{marker}"] = parent.astype(int)
                    table[marker] = author.astype(int)
                    table[f"projected_{marker}"] = projected.astype(int)
                    references.append(dict(sample=sid, population=population, marker=marker,
                        importedCount=int(author.sum()), projectedCount=int(projected.sum()),
                        rawRectangleCount=int(raw_rectangle.sum()),
                        rawRectangleMismatch=int((raw_rectangle != author).sum()),
                        rawParentMismatch=int((raw_parent != parent).sum()),
                        projectedOutsideAuthor=int((projected & ~author).sum()),
                        cytokineLowerRaw=lower, cytokineUpperRaw=upper,
                        cd44LowerRaw=limits[CD44]["min"], cd44UpperRaw=limits[CD44]["max"],
                        aboveAuthorUpperCount=int((parent & (cytokine >= upper)).sum()),
                        importedParent=int(parent.sum()),
                        lineageCd44Parent=int((lineage & in_interval).sum()),
                        importedLineage=int(lineage.sum()), flowjoLineage=stored[(path[-1], path[:-1])],
                        importedNonNaive=int(non_naive.sum()), flowjoNonNaive=stored[("Non-naive", path)],
                        flowjoCount=stored[(GATES[marker], child_path)]))
                table.to_csv(output_dir / filename, index=False)
                parents.append(dict(sample=sid, population=population, nCell=int(mask.sum()),
                                    nNonNaive=int(non_naive[mask].sum()),
                                    file=filename))
            print(f"Audited {sid}", flush=True)
    pd.DataFrame(audits).to_csv(output_dir / "gate-validation.csv", index=False)
    pd.DataFrame(parents).to_csv(output_dir / "populations.csv", index=False)
    pd.DataFrame(references).to_csv(output_dir / "references.csv", index=False)
    (output_dir / "importer.json").write_text(json.dumps({"FlowKit": fk.__version__}) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("raw_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    export(args.raw_dir, args.output_dir)
