"""Reconstruct OMIP-111 parent and cytokine gate membership using FlowKit 1.3.2.

Export original FCS event indices, never the clipped FlowJo display values.
Audit every stored gate count. The R analysis applies a common raw-data transform.
"""
import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import flowkit as fk
import numpy as np
import pandas as pd

CHANNELS = {"IFNg": "AF488-A", "IL2": "PE-A", "TNF": "PE-Cy7-A", "IL4_5": "BV421-A", "IL17A": "R718-A"}
GATES = {"IFNg": "IFNg+", "IL2": "IL-2", "TNF": "TNF+",
         "IL4_5": "Th2 Cells", "IL17A": "Th17 Cells"}


def stored_gates(node, path=("root",), definitions=None):
    result = {}
    for child in node.findall("./Subpopulations/Population"):
        name = child.attrib["name"]
        result[(name, path)] = int(child.attrib["count"])
        if definitions is not None:
            definitions[(name, path)] = child
        result.update(stored_gates(child, path + (name,), definitions))
    return result


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
            for name, path in wsp.get_gate_ids(sid):
                observed = int(wsp.get_gate_membership(sid, name, gate_path=path).sum())
                expected = stored[(name, path)]
                audits.append(dict(sample=sid, pop="/" + "/".join((path + (name,))[1:]),
                                   count_import=observed, count_flowjo=expected, difference=observed-expected))
            for population in ("CD4", "CD8"):
                paths = [p for p in wsp.find_matching_gate_paths(sid, "Non-naive")
                         if p[-1] == f"{population}+ T Cells"]
                if len(paths) != 1:
                    raise ValueError("Ambiguous non-naive parent")
                path = paths[0]
                mask = wsp.get_gate_membership(sid, "Non-naive", gate_path=path)
                markers = list(GATES) if population == "CD4" else list(GATES)[:3]
                table = pd.DataFrame({"eventIndex": np.flatnonzero(mask) + 1})
                filename = sid.replace(".fcs", "") + f"-{population}.csv"
                parent_flowjo = stored[("Non-naive", path)]
                for marker in markers:
                    child_path = path + ("Non-naive",)
                    membership = wsp.get_gate_membership(sid, GATES[marker], gate_path=child_path)[mask]
                    definition = definitions[(GATES[marker], child_path)]
                    dimensions = [d for d in definition.iter() if d.tag.rsplit("}", 1)[-1] == "dimension"]
                    cutoffs = {}
                    for dimension in dimensions:
                        channel_nodes = [n for n in dimension.iter() if n.tag.rsplit("}", 1)[-1] == "fcs-dimension"]
                        if len(channel_nodes) != 1:
                            raise ValueError("Unexpected author gate dimension")
                        channel = next(v for k, v in channel_nodes[0].attrib.items() if k.rsplit("}", 1)[-1] == "name")
                        cutoffs[channel] = {k.rsplit("}", 1)[-1]: float(v) for k, v in dimension.attrib.items()}
                    if set(cutoffs) != {CHANNELS[marker], "[RealBlue 744]-A"}:
                        raise ValueError("Author cytokine gate is not cytokine-by-CD44")
                    cytokine = raw_events[mask, sample.get_channel_index(CHANNELS[marker])]
                    cutoff = cutoffs[CHANNELS[marker]]["min"]
                    upper = cutoffs[CHANNELS[marker]]["max"]
                    projected = cytokine > cutoff
                    table[marker] = projected.astype(int)
                    table[f"author_{marker}"] = membership.astype(int)
                    references.append(dict(sample=sid, population=population, marker=marker,
                        importedCount=int(membership.sum()), projectedCount=int(projected.sum()),
                        cytokineLowerRaw=cutoff, cytokineUpperRaw=upper,
                        aboveAuthorUpperCount=int((cytokine > upper).sum()),
                        importedParent=int(mask.sum()),
                        flowjoCount=stored[(GATES[marker], child_path)], flowjoParent=parent_flowjo))
                table.to_csv(output_dir / filename, index=False)
                parents.append(dict(sample=sid, population=population, nCell=int(mask.sum()), file=filename))
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
