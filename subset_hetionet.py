"""
Hetionet Subset Builder v2
==========================
Builds a disease-focused subset with ALL 11 entity types represented.

Usage:
    python subset_hetionet.py

Input:  data/hetionet-v1.0.json.bz2
Output: data/hetionet_subset.json
"""

import json
import bz2
import os
from collections import defaultdict

# ── 1. Load the full graph ────────────────────────────────────────
input_path = os.path.join("data", "hetionet-v1.0.json.bz2")
print(f"Loading {input_path} (this may take 30-60 seconds)...")
with bz2.open(input_path, "rt") as f:
    data = json.load(f)

print(f"Full graph: {len(data['nodes'])} nodes, {len(data['edges'])} edges")

# ── 2. Index nodes ────────────────────────────────────────────────
nodes_by_id = {}
for node in data["nodes"]:
    nid = f"{node['kind']}::{node['identifier']}"
    nodes_by_id[nid] = node

kind_counts = defaultdict(int)
for node in data["nodes"]:
    kind_counts[node["kind"]] += 1
print("\nFull graph node types:")
for kind, count in sorted(kind_counts.items(), key=lambda x: -x[1]):
    print(f"  {kind}: {count}")

# ── 3. Build adjacency ───────────────────────────────────────────
adj = defaultdict(set)
edges_by_node = defaultdict(list)

for edge in data["edges"]:
    src = f"{edge['source_id'][0]}::{edge['source_id'][1]}"
    tgt = f"{edge['target_id'][0]}::{edge['target_id'][1]}"
    adj[src].add(tgt)
    adj[tgt].add(src)
    edges_by_node[src].append(edge)
    edges_by_node[tgt].append(edge)

# ── 4. Select seed diseases ───────────────────────────────────────
diseases_with_treatments = set()
for edge in data["edges"]:
    if edge["kind"] == "treats":
        disease_id = f"{edge['target_id'][0]}::{edge['target_id'][1]}"
        diseases_with_treatments.add(disease_id)

print(f"\nDiseases with treatments: {len(diseases_with_treatments)}")

seed_keywords = [
    "alzheimer", "parkinson", "epilepsy", "schizophrenia",
    "breast cancer", "lung cancer", "leukemia",
    "diabetes", "obesity", "hypertension",
    "asthma", "rheumatoid arthritis", "depression",
]

seed_diseases = []
for nid, node in nodes_by_id.items():
    if node["kind"] == "Disease" and nid in diseases_with_treatments:
        name_lower = node["name"].lower()
        for kw in seed_keywords:
            if kw in name_lower:
                seed_diseases.append(nid)
                break

if len(seed_diseases) < 10:
    disease_degree = {}
    for nid in diseases_with_treatments:
        disease_degree[nid] = len(adj.get(nid, set()))
    top_diseases = sorted(disease_degree.items(), key=lambda x: -x[1])[:15]
    for nid, _ in top_diseases:
        if nid not in seed_diseases:
            seed_diseases.append(nid)

seed_diseases = seed_diseases[:15]
seed_ids = set(seed_diseases)

print(f"\nSeed diseases ({len(seed_diseases)}):")
for sid in seed_diseases:
    print(f"  {nodes_by_id[sid]['name']} (degree: {len(adj.get(sid, set()))})")

# ── 5. BFS hop 1: all direct neighbors of seeds ──────────────────
subset_nodes = set(seed_diseases)
hop1 = set()
for node in seed_diseases:
    for neighbor in adj.get(node, set()):
        hop1.add(neighbor)
subset_nodes.update(hop1)
print(f"\nAfter hop 1: {len(subset_nodes)} nodes")

# ── 6. BFS hop 2: expand from Compound and Gene nodes ────────────
# Now we ALLOW all entity types in hop 2, but cap each type
hop2_all = defaultdict(set)  # kind -> set of node ids

for node in hop1:
    node_kind = nodes_by_id.get(node, {}).get("kind", "")
    if node_kind in ("Compound", "Gene"):
        for neighbor in adj.get(node, set()):
            if neighbor in subset_nodes:
                continue
            nkind = nodes_by_id.get(neighbor, {}).get("kind", "")
            hop2_all[nkind].add(neighbor)

print("\nHop 2 candidates by type (before capping):")
for kind, nids in sorted(hop2_all.items(), key=lambda x: -len(x[1])):
    print(f"  {kind}: {len(nids)}")

# Cap per type: keep the most connected nodes of each type
hop2_caps = {
    "Disease": 999,              # keep all (small)
    "Compound": 200,             # drugs
    "Side Effect": 200,          # clinical
    "Symptom": 100,              # clinical
    "Pharmacologic Class": 999,  # keep all (small)
    "Gene": 200,                 # already have many from hop 1
    "Anatomy": 100,              # body locations
    "Biological Process": 60,    # key mechanisms
    "Molecular Function": 40,    # molecular roles
    "Cellular Component": 30,    # cell structures
    "Pathway": 80,               # biological pathways
}

for kind, nids in hop2_all.items():
    cap = hop2_caps.get(kind, 50)
    if len(nids) <= cap:
        subset_nodes.update(nids)
    else:
        # Keep most-connected within the existing subset
        scored = []
        for nid in nids:
            conn_in_subset = sum(1 for n in adj.get(nid, set()) if n in subset_nodes)
            scored.append((conn_in_subset, nid))
        scored.sort(key=lambda x: -x[0])
        for _, nid in scored[:cap]:
            subset_nodes.add(nid)

print(f"\nAfter hop 2 (capped): {len(subset_nodes)} nodes")

# ── 7. Overall cap ────────────────────────────────────────────────
MAX_NODES = 3000
if len(subset_nodes) > MAX_NODES:
    node_degree = {}
    for nid in subset_nodes:
        node_degree[nid] = sum(1 for n in adj.get(nid, set()) if n in subset_nodes)

    keep = set(seed_diseases)
    remaining = [(nid, node_degree[nid]) for nid in subset_nodes if nid not in keep]
    remaining.sort(key=lambda x: -x[1])
    for nid, _ in remaining:
        if len(keep) >= MAX_NODES:
            break
        keep.add(nid)
    subset_nodes = keep
    print(f"After overall pruning: {len(subset_nodes)} nodes")

# ── 8. Collect edges ──────────────────────────────────────────────
subset_edges = []
seen_edges = set()

for edge in data["edges"]:
    src = f"{edge['source_id'][0]}::{edge['source_id'][1]}"
    tgt = f"{edge['target_id'][0]}::{edge['target_id'][1]}"
    if src in subset_nodes and tgt in subset_nodes:
        edge_key = (src, tgt, edge["kind"])
        if edge_key not in seen_edges:
            seen_edges.add(edge_key)
            subset_edges.append(edge)

print(f"Subset edges (unfiltered): {len(subset_edges)}")

# ── 9. Edge budget to keep D3 performant ──────────────────────────
edge_budgets = {
    "treats": 9999,
    "palliates": 9999,
    "presents": 9999,       # disease-symptom (important!)
    "includes": 9999,       # pharmacologic class (small)
    "localizes": 9999,
    "binds": 800,
    "associates": 800,
    "causes": 600,
    "resembles": 150,
    "participates": 500,    # gene-pathway/BP/MF/CC (increased!)
    "expresses": 300,
    "upregulates": 200,
    "downregulates": 200,
    "regulates": 150,
    "interacts": 150,
    "covaries": 80,
}

# Score edges: seed-connected first, then by node degree
node_degree = defaultdict(int)
for e in subset_edges:
    src = f"{e['source_id'][0]}::{e['source_id'][1]}"
    tgt = f"{e['target_id'][0]}::{e['target_id'][1]}"
    node_degree[src] += 1
    node_degree[tgt] += 1

by_type = defaultdict(list)
for e in subset_edges:
    by_type[e["kind"]].append(e)

filtered_edges = []
for kind, elist in by_type.items():
    budget = edge_budgets.get(kind, 80)
    scored = []
    for e in elist:
        src = f"{e['source_id'][0]}::{e['source_id'][1]}"
        tgt = f"{e['target_id'][0]}::{e['target_id'][1]}"
        score = 0
        if src in seed_ids or tgt in seed_ids:
            score += 10000
        score += node_degree[src] + node_degree[tgt]
        scored.append((score, e))
    scored.sort(key=lambda x: -x[0])
    filtered_edges.extend([e for _, e in scored[:budget]])

print(f"Subset edges (after budget): {len(filtered_edges)}")

# ── 10. Prune disconnected nodes ─────────────────────────────────
connected = set()
for e in filtered_edges:
    connected.add(f"{e['source_id'][0]}::{e['source_id'][1]}")
    connected.add(f"{e['target_id'][0]}::{e['target_id'][1]}")
for sid in seed_ids:
    connected.add(sid)

filtered_nodes = [nid for nid in subset_nodes if nid in connected]

# ── 11. Build output ─────────────────────────────────────────────
layer_map = {
    "Disease": 1, "Symptom": 1, "Side Effect": 1,
    "Compound": 2, "Pharmacologic Class": 2,
    "Gene": 3, "Anatomy": 3,
    "Biological Process": 4, "Pathway": 4,
    "Molecular Function": 4, "Cellular Component": 4,
}
layer_names = {1: "Clinical", 2: "Pharmaceutical", 3: "Genomic", 4: "Biological"}

output_nodes = []
for nid in filtered_nodes:
    node = nodes_by_id[nid]
    output_nodes.append({
        "id": nid,
        "name": node["name"],
        "kind": node["kind"],
        "layer": layer_map.get(node["kind"], 4),
        "layer_name": layer_names.get(layer_map.get(node["kind"], 4), "Other"),
        "identifier": str(node["identifier"]),
    })

output_edges = []
for edge in filtered_edges:
    output_edges.append({
        "source": f"{edge['source_id'][0]}::{edge['source_id'][1]}",
        "target": f"{edge['target_id'][0]}::{edge['target_id'][1]}",
        "kind": edge["kind"],
        "direction": edge["direction"],
    })

output = {
    "nodes": output_nodes,
    "edges": output_edges,
    "seed_diseases": [nodes_by_id[sid]["name"] for sid in seed_diseases],
    "meta": {
        "total_nodes": len(output_nodes),
        "total_edges": len(output_edges),
        "layers": layer_names,
    },
}

# ── 12. Summary ───────────────────────────────────────────────────
print(f"\n{'='*55}")
print(f"  SUBSET SUMMARY")
print(f"{'='*55}")
print(f"  Nodes: {len(output_nodes)}")
print(f"  Edges: {len(output_edges)}")

# Check all 11 types
all_kinds = set()
layer_counts = defaultdict(lambda: defaultdict(int))
for n in output_nodes:
    layer_counts[n["layer_name"]][n["kind"]] += 1
    all_kinds.add(n["kind"])

print(f"\n  Nodes by layer:")
for layer_name in ["Clinical", "Pharmaceutical", "Genomic", "Biological"]:
    kinds = layer_counts.get(layer_name, {})
    total = sum(kinds.values())
    detail = ", ".join(f"{k}: {v}" for k, v in sorted(kinds.items()))
    print(f"    {layer_name}: {total} ({detail})")

expected_kinds = {"Disease", "Symptom", "Side Effect", "Compound",
                  "Pharmacologic Class", "Gene", "Anatomy",
                  "Biological Process", "Pathway", "Molecular Function",
                  "Cellular Component"}
missing = expected_kinds - all_kinds
if missing:
    print(f"\n  ⚠ MISSING entity types: {missing}")
else:
    print(f"\n  ✓ All 11 entity types represented!")

print(f"\n  Edges by type:")
edge_type_counts = defaultdict(int)
for e in output_edges:
    edge_type_counts[e["kind"]] += 1
for kind, count in sorted(edge_type_counts.items(), key=lambda x: -x[1]):
    print(f"    {kind}: {count}")

# ── 13. Save ──────────────────────────────────────────────────────
output_path = os.path.join("data", "hetionet_subset.json")
with open(output_path, "w") as f:
    json.dump(output, f, separators=(",", ":"))

file_size = os.path.getsize(output_path) / 1024 / 1024
print(f"\n  Saved to {output_path} ({file_size:.1f} MB)")
print(f"  Upload this file to Claude to continue.")