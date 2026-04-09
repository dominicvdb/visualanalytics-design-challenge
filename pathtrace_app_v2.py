import marimo

__generated_with = "0.21.1"
app = marimo.App(width="full", app_title="Biomedical Knowledge Graph Explorer")


@app.cell
def _():
    import json as json_lib
    import marimo as mo
    import networkx as nx
    import traitlets
    from collections import Counter, defaultdict

    try:
        import anywidget
    except ImportError:
        anywidget = None

    return Counter, defaultdict, json_lib, mo, nx, traitlets, anywidget


@app.cell
def _(json_lib):
    with open("data/hetionet_subset.json", "r", encoding="utf-8") as _f:
        graph_data = json_lib.load(_f)

    print(f"Loaded: {graph_data['meta']['total_nodes']} nodes, {graph_data['meta']['total_edges']} edges")
    print(f"Seed diseases: {', '.join(graph_data['seed_diseases'][:8])}")
    return (graph_data,)


@app.cell
def _(Counter, defaultdict, graph_data, nx):
    G = nx.Graph()
    nodes_by_id = {}
    edge_types_by_node = defaultdict(set)
    neighbor_layers = defaultdict(set)
    cross_layer_counts = Counter()

    for _node in graph_data["nodes"]:
        G.add_node(_node["id"], **_node)
        nodes_by_id[_node["id"]] = dict(_node)

    seed_name_set = set(graph_data["seed_diseases"])
    seed_ids = [_node["id"] for _node in graph_data["nodes"] if _node["name"] in seed_name_set]

    for _edge in graph_data["edges"]:
        _src = _edge["source"]
        _tgt = _edge["target"]
        G.add_edge(_src, _tgt, kind=_edge["kind"], direction=_edge["direction"])
        edge_types_by_node[_src].add(_edge["kind"])
        edge_types_by_node[_tgt].add(_edge["kind"])

        _src_layer = nodes_by_id[_src]["layer"]
        _tgt_layer = nodes_by_id[_tgt]["layer"]
        neighbor_layers[_src].add(_tgt_layer)
        neighbor_layers[_tgt].add(_src_layer)
        if _src_layer != _tgt_layer:
            cross_layer_counts[_src] += 1
            cross_layer_counts[_tgt] += 1

    try:
        seed_distances = nx.multi_source_dijkstra_path_length(G, seed_ids) if seed_ids else {}
    except Exception:
        seed_distances = {}

    for _nid, _node in nodes_by_id.items():
        _degree = G.degree(_nid)
        _unique_relations = len(edge_types_by_node[_nid])
        _unique_layers = len(neighbor_layers[_nid])
        _seed_distance = seed_distances.get(_nid)
        _seed_bonus = 0.0
        if _seed_distance is not None:
            _seed_bonus = max(0, 4 - min(_seed_distance, 4))

        _discovery_score = (
            _unique_layers * 3.0
            + _unique_relations * 2.0
            + cross_layer_counts[_nid] * 1.5
            + min(_degree, 20) * 0.4
            + _seed_bonus
        )

        _node["degree"] = _degree
        _node["edge_type_count"] = _unique_relations
        _node["neighbor_layer_count"] = _unique_layers
        _node["cross_layer_count"] = cross_layer_counts[_nid]
        _node["seed_distance"] = _seed_distance
        _node["discovery_score"] = round(_discovery_score, 2)
        _node["is_seed"] = _node["name"] in seed_name_set

    _component_sizes = sorted((len(_c) for _c in nx.connected_components(G)), reverse=True)
    _largest_component = _component_sizes[0] if _component_sizes else 0
    _degree_zero = sum(1 for _nid in nodes_by_id if nodes_by_id[_nid]["degree"] == 0)
    _degree_one = sum(1 for _nid in nodes_by_id if nodes_by_id[_nid]["degree"] == 1)

    global_metrics = {
        "largest_component": _largest_component,
        "largest_component_pct": round((_largest_component / max(G.number_of_nodes(), 1)) * 100, 1),
        "degree_zero": _degree_zero,
        "degree_one": _degree_one,
        "seed_count": len(seed_ids),
    }

    print(f"Graph built: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    print(
        f"Largest connected component: {global_metrics['largest_component']} "
        f"({global_metrics['largest_component_pct']}%)"
    )
    return G, global_metrics, nodes_by_id, seed_ids


@app.cell
def _(mo):
    search_input = mo.ui.text(
        placeholder="Search node (e.g. Metformin, Alzheimer's, BRCA1)...",
        label="Search",
        full_width=True,
    )
    return (search_input,)


@app.cell
def _(mo):
    edge_preset = mo.ui.dropdown(
        options=[
            "Cross-layer discovery",
            "Clinical focus",
            "Biological mechanisms",
            "All relations",
        ],
        value="Cross-layer discovery",
        label="Relation preset",
    )
    return (edge_preset,)


@app.cell
def _(mo):
    sort_mode = mo.ui.dropdown(
        options=[
            "Discovery score",
            "Degree",
            "Alphabetical",
            "Seed proximity",
        ],
        value="Discovery score",
        label="Sort within layers",
    )
    return (sort_mode,)



@app.cell
def _(edge_preset, graph_data, global_metrics, json_lib, nodes_by_id, search_input, sort_mode):
    layer_config = {
        1: {"label": "Clinical", "color": "#E74C3C", "sublabel": "Disease · Symptom · Side Effect"},
        2: {"label": "Pharmaceutical", "color": "#3498DB", "sublabel": "Compound · Pharmacologic Class"},
        3: {"label": "Genomic", "color": "#2ECC71", "sublabel": "Gene · Anatomy"},
        4: {"label": "Biological", "color": "#9B59B6", "sublabel": "Pathway · Process · Function · Component"},
    }

    kind_colors = {
        "Disease": "#E74C3C", "Symptom": "#E67E22", "Side Effect": "#F39C12",
        "Compound": "#3498DB", "Pharmacologic Class": "#2980B9",
        "Gene": "#2ECC71", "Anatomy": "#27AE60",
        "Biological Process": "#9B59B6", "Pathway": "#8E44AD",
        "Molecular Function": "#A569BD", "Cellular Component": "#7D3C98",
    }

    edge_colors = {
        "treats": "#27AE60", "palliates": "#85C1E9",
        "binds": "#F4D03F", "associates": "#E74C3C",
        "causes": "#E74C3C", "resembles": "#95A5A6",
        "expresses": "#82E0AA", "participates": "#BB8FCE",
        "upregulates": "#F5B041", "downregulates": "#5DADE2",
        "regulates": "#EB984E", "interacts": "#AEB6BF",
        "localizes": "#73C6B6", "covaries": "#D7BDE2",
    }

    preset_map = {
        "Cross-layer discovery": ["treats", "binds", "associates", "participates", "localizes", "expresses"],
        "Clinical focus": ["treats", "palliates", "associates", "causes", "resembles"],
        "Biological mechanisms": ["binds", "participates", "localizes", "expresses", "upregulates", "downregulates", "regulates", "interacts"],
        "All relations": list(edge_colors.keys()),
    }

    active_types = [
        _edge_type for _edge_type in preset_map.get(edge_preset.value, preset_map["Cross-layer discovery"])
        if _edge_type
    ]
    active_type_set = set(active_types)

    layer_groups = {1: [], 2: [], 3: [], 4: []}
    for _node in graph_data["nodes"]:
        layer_groups[_node["layer"]].append(_node)

    def _sort_key(_node):
        _metrics = nodes_by_id[_node["id"]]
        if sort_mode.value == "Degree":
            return (-_metrics["degree"], _node["name"].lower())
        if sort_mode.value == "Alphabetical":
            return (_node["name"].lower(),)
        if sort_mode.value == "Seed proximity":
            _dist = _metrics["seed_distance"] if _metrics["seed_distance"] is not None else 999
            return (_dist, -_metrics["discovery_score"], _node["name"].lower())
        return (-_metrics["discovery_score"], -_metrics["degree"], _node["name"].lower())

    plot_nodes = []
    discovery_leads = []
    for _layer_id, _nlist in layer_groups.items():
        _sorted_nodes = sorted(_nlist, key=_sort_key)
        for _idx, _node in enumerate(_sorted_nodes):
            _metrics = nodes_by_id[_node["id"]]
            _payload = {
                "id": _node["id"],
                "name": _node["name"],
                "kind": _node["kind"],
                "layer": _node["layer"],
                "layerName": layer_config[_node["layer"]]["label"],
                "color": kind_colors.get(_node["kind"], "#999"),
                "degree": _metrics["degree"],
                "rank": _idx,
                "layerSize": len(_sorted_nodes),
                "isSeed": _metrics["is_seed"],
                "discoveryScore": _metrics["discovery_score"],
                "edgeTypeCount": _metrics["edge_type_count"],
                "neighborLayerCount": _metrics["neighbor_layer_count"],
                "crossLayerCount": _metrics["cross_layer_count"],
                "seedDistance": _metrics["seed_distance"],
            }
            plot_nodes.append(_payload)
            discovery_leads.append(_payload)

    discovery_leads = sorted(
        discovery_leads,
        key=lambda _n: (-_n["discoveryScore"], -_n["crossLayerCount"], -_n["degree"]),
    )[:8]

    plot_edges = []
    for _edge in graph_data["edges"]:
        if _edge["kind"] not in active_type_set:
            continue
        plot_edges.append({
            "source": _edge["source"],
            "target": _edge["target"],
            "kind": _edge["kind"],
            "direction": _edge["direction"],
            "color": edge_colors.get(_edge["kind"], "#BDBDBD"),
        })

    payload = {
        "nodes": plot_nodes,
        "edges": plot_edges,
        "layerConfig": layer_config,
        "edgeColors": edge_colors,
        "search": (search_input.value or "").strip(),
        "discoveryLeads": discovery_leads,
        "globalMetrics": global_metrics,
        "activeEdgeCount": len(plot_edges),
        "activeTypes": active_types,
        "viewLabel": edge_preset.value,
        "sortLabel": sort_mode.value,
    }

    payload_json = json_lib.dumps(payload)
    return active_types, payload_json


@app.cell
def _(active_types, edge_preset, global_metrics, mo, sort_mode):
    summary_md = mo.md(
        f"""
# Biomedical Knowledge Graph Explorer

Exploring **4 biomedical layers** with **{global_metrics['seed_count']}** seed diseases.  
View: **{edge_preset.value}** · Sort: **{sort_mode.value}** · Showing **{len(active_types)}** active relation types.

Connectedness: **{global_metrics['largest_component_pct']}%** in largest component · degree-0 nodes: **{global_metrics['degree_zero']}**.
"""
    )
    return (summary_md,)


@app.cell
def _(anywidget, mo):
    anywidget_install_help = None
    if anywidget is None:
        anywidget_install_help = mo.md(
            """
**D3 widget dependency missing.** Install **anywidget** in your environment and rerun:

```bash
pip install anywidget
```

This app uses a D3-based anywidget instead of `mo.iframe(...)` so the JavaScript runs reliably inside marimo.
"""
        )
    return (anywidget_install_help,)


@app.cell
def _(anywidget, traitlets):
    if anywidget is None:
        GraphWidget = None
    else:
        class GraphWidget(anywidget.AnyWidget):
            payload = traitlets.Unicode("").tag(sync=True)
            _css = """
            .pt-root { font-family: Inter, system-ui, sans-serif; color: #ddd; }
            .pt-top-bar { display: flex; gap: 12px; margin-bottom: 8px; align-items: flex-start; }
            .pt-top-leads { flex: 1; display: flex; gap: 8px; flex-wrap: wrap; align-items: flex-start; }
            .pt-lead-card { border: 1px solid #263247; background: rgba(10,13,20,0.92); border-radius: 10px; padding: 8px 12px; cursor: pointer; min-width: 160px; max-width: 200px; flex: 1; }
            .pt-lead-card:hover { background: rgba(138,164,255,0.10); border-color: #3a5070; }
            .pt-lead-card .pt-lead-name { font-weight: 700; font-size: 12px; }
            .pt-lead-card .pt-lead-detail { font-size: 10px; color: #7d8596; margin-top: 2px; }
            .pt-top-label { color: #8aa4ff; font-weight: 700; font-size: 10px; text-transform: uppercase; letter-spacing: 0.08em; padding-top: 10px; white-space: nowrap; }
            .pt-shell { display: flex; width: 100%; height: 860px; background: #070a11; border: 1px solid #1e2633; border-radius: 12px; overflow: hidden; }
            .pt-graph { flex: 1; position: relative; overflow: hidden; background: #070a11; }
            .pt-side { width: 340px; background: rgba(10,13,20,0.98); border-left: 1px solid #1e2633; overflow-y: auto; }
            .pt-header { padding: 14px 16px 10px; border-bottom: 1px solid #1e2633; position: sticky; top: 0; background: rgba(10,13,20,0.98); z-index: 2; }
            .pt-header h3 { margin: 0; font-size: 16px; }
            .pt-sub { color: #7d8596; font-size: 11px; margin-top: 3px; }
            .pt-body { padding: 14px 16px 24px; }
            .pt-section { margin-bottom: 14px; }
            .pt-section-title { color: #8aa4ff; font-weight: 700; font-size: 11px; text-transform: uppercase; letter-spacing: 0.08em; margin-bottom: 6px; }
            .pt-card { border: 1px solid #263247; background: rgba(255,255,255,0.02); border-radius: 10px; padding: 10px 12px; margin-bottom: 8px; }
            .pt-chip-wrap { display: flex; flex-wrap: wrap; gap: 6px; }
            .pt-chip { display: inline-block; padding: 4px 8px; border-radius: 999px; border: 1px solid #34445f; background: rgba(255,255,255,0.03); color: #d8deea; font-size: 11px; cursor: pointer; }
            .pt-chip:hover { background: rgba(138,164,255,0.12); }
            .pt-overlay { position: absolute; background: rgba(10,13,20,0.92); border: 1px solid #1f2b3d; border-radius: 10px; padding: 8px 10px; font-size: 11px; }
            .pt-stats { top: 12px; left: 12px; color: #b8c2d4; }
            .pt-legend { top: 12px; right: 12px; width: 220px; max-height: 300px; overflow-y: auto; }
            .pt-overlay-title { color: #8aa4ff; font-weight: 700; text-transform: uppercase; letter-spacing: 0.08em; font-size: 10px; margin-bottom: 6px; }
            .pt-legend-row { display: flex; align-items: center; gap: 6px; margin: 3px 0; color: #b8c2d4; }
            .pt-swatch { width: 10px; height: 10px; border-radius: 2px; flex-shrink: 0; }
            .pt-tooltip { position: fixed; pointer-events: none; z-index: 9999; display: none; background: rgba(10,13,20,0.96); border: 1px solid #34445f; color: #e2e8f0; padding: 8px 10px; border-radius: 8px; font-size: 12px; max-width: 280px; box-shadow: 0 10px 30px rgba(0,0,0,0.35); }
            .pt-empty { color: #7d8596; }
            .pt-stat-grid { display: grid; grid-template-columns: repeat(3, 1fr); gap: 8px; margin-top: 8px; }
            .pt-mini { text-align: center; background: rgba(255,255,255,0.03); border: 1px solid #263247; border-radius: 8px; padding: 6px; }
            .pt-mini .v { font-size: 15px; font-weight: 700; color: #f5f7fb; }
            .pt-mini .l { font-size: 9px; color: #7d8596; text-transform: uppercase; letter-spacing: 0.06em; }
            svg { width: 100%; height: 100%; display: block; }
            """
            _esm = r'''
            import * as d3 from "https://esm.sh/d3@7";

            function render({ model, el }) {
              const draw = () => {
                const state = JSON.parse(model.get("payload") || "{}");
                const nodes = state.nodes || [];
                const edges = state.edges || [];
                const layerConfig = state.layerConfig || {};
                const searchTerm = (state.search || "").toLowerCase();
                const leads = state.discoveryLeads || [];
                const root = document.createElement("div");
                root.className = "pt-root";
                root.innerHTML = `
                  <div class="pt-top-bar">
                    <div class="pt-top-label">Top discovery leads</div>
                    <div class="pt-top-leads" id="pt-leads-bar"></div>
                  </div>
                  <div class="pt-shell">
                    <div class="pt-graph">
                      <div class="pt-overlay pt-stats"></div>
                      <div class="pt-overlay pt-legend"></div>
                    </div>
                    <div class="pt-side">
                      <div class="pt-header">
                        <h3>Discovery panel</h3>
                        <div class="pt-sub">Click a node to inspect neighbors, bridge paths, and relation types.</div>
                      </div>
                      <div class="pt-body"><div class="pt-empty">Select a node to begin.</div></div>
                    </div>
                  </div>
                  <div class="pt-tooltip"></div>
                `;
                el.replaceChildren(root);

                const graphArea = root.querySelector(".pt-graph");
                const sideBody = root.querySelector(".pt-body");
                const statsBox = root.querySelector(".pt-stats");
                const legendBox = root.querySelector(".pt-legend");
                const leadsBox = root.querySelector("#pt-leads-bar");
                const tooltip = root.querySelector(".pt-tooltip");

                const width = Math.max(graphArea.clientWidth || 1000, 900);
                const height = 920;
                const margin = { top: 60, right: 36, bottom: 30, left: 60 };
                const plotW = width - margin.left - margin.right;
                const plotH = height - margin.top - margin.bottom;
                const layerH = plotH / 4;

                const nodeMap = new Map(nodes.map((n) => [n.id, { ...n }]));
                const adj = new Map();
                for (const _edge of edges) {
                  if (!adj.has(_edge.source)) adj.set(_edge.source, []);
                  if (!adj.has(_edge.target)) adj.set(_edge.target, []);
                  adj.get(_edge.source).push(_edge);
                  adj.get(_edge.target).push(_edge);
                }

                // ── Per-layer force simulation for overlap-free layout ──
                // Based on constrained force-directed placement (Dwyer et al., 2006)
                // with d3.forceCollide for node overlap removal.
                const layerGroups = {1:[], 2:[], 3:[], 4:[]};
                for (const _node of nodes) {
                  if (layerGroups[_node.layer]) layerGroups[_node.layer].push(_node);
                }

                function nodeRadius(n) {
                  return n.isSeed ? 7 : Math.max(2.5, Math.min(5.5, Math.sqrt(n.degree) * 0.8));
                }

                for (const lid of [1,2,3,4]) {
                  const group = layerGroups[lid];
                  if (!group.length) continue;

                  const bandY0 = margin.top + (lid - 1) * layerH + 40;
                  const bandY1 = margin.top + lid * layerH - 18;
                  const bandH = bandY1 - bandY0;
                  const bandCY = (bandY0 + bandY1) / 2;

                  // Initial positions: sorted by Python-computed rank (respects sort mode)
                  group.sort((a, b) => a.rank - b.rank);
                  const cols = Math.max(1, Math.ceil(Math.sqrt(group.length * (plotW / bandH))));
                  const rows = Math.max(1, Math.ceil(group.length / cols));
                  const cellW = plotW / cols;
                  const cellH = bandH / rows;

                  group.forEach((_node, i) => {
                    const col = i % cols;
                    const row = Math.floor(i / cols);
                    _node.x = margin.left + col * cellW + cellW / 2;
                    _node.y = bandY0 + row * cellH + cellH / 2;
                    _node._r = nodeRadius(_node) + 1.5;
                  });

                  // Run constrained force simulation
                  const sim = d3.forceSimulation(group)
                    .force("collide", d3.forceCollide(d => d._r + 1).strength(0.8).iterations(3))
                    .force("x", d3.forceX(d => d.x).strength(0.3))
                    .force("y", d3.forceY(d => d.y).strength(0.3))
                    .stop();

                  // Run ticks synchronously
                  const ticks = Math.min(80, 30 + group.length / 10);
                  for (let t = 0; t < ticks; t++) sim.tick();

                  // Clamp to layer bounds
                  for (const _node of group) {
                    _node.x = Math.max(margin.left + 5, Math.min(margin.left + plotW - 5, _node.x));
                    _node.y = Math.max(bandY0 + 3, Math.min(bandY1 - 3, _node.y));
                    nodeMap.set(_node.id, _node);
                  }
                }

                const svg = d3.select(graphArea)
                  .append("svg")
                  .attr("viewBox", [0, 0, width, height])
                  .attr("preserveAspectRatio", "xMidYMid meet");

                const zoomLayer = svg.append("g");
                svg.call(
                  d3.zoom()
                    .scaleExtent([0.35, 7])
                    .on("zoom", (event) => {
                      zoomLayer.attr("transform", event.transform);
                    })
                );

                for (const lid of [1, 2, 3, 4]) {
                  const cfg = layerConfig[lid] || { label: `Layer ${lid}`, color: "#888", sublabel: "" };
                  const y0 = margin.top + (lid - 1) * layerH;
                  zoomLayer.append("rect")
                    .attr("x", margin.left - 20)
                    .attr("y", y0)
                    .attr("width", plotW + 40)
                    .attr("height", layerH - 10)
                    .attr("rx", 8)
                    .attr("fill", "rgba(255,255,255,0.02)")
                    .attr("stroke", "rgba(255,255,255,0.06)");
                  zoomLayer.append("text")
                    .attr("x", margin.left - 10)
                    .attr("y", y0 + 18)
                    .attr("fill", cfg.color)
                    .attr("font-size", 13)
                    .attr("font-weight", 700)
                    .text(`Layer ${lid}: ${cfg.label}`);
                  zoomLayer.append("text")
                    .attr("x", margin.left - 10)
                    .attr("y", y0 + 33)
                    .attr("fill", "#6b7280")
                    .attr("font-size", 10)
                    .text(cfg.sublabel || "");
                }

                const edgeGroup = zoomLayer.append("g");
                const nodeGroup = zoomLayer.append("g");

                const edgeSel = edgeGroup.selectAll("path")
                  .data(edges)
                  .enter()
                  .append("path")
                  .attr("fill", "none")
                  .attr("stroke", (d) => d.color || "#bdbdbd")
                  .attr("stroke-width", 0.6)
                  .attr("opacity", 0.10)
                  .attr("d", (d) => {
                    const s = nodeMap.get(d.source);
                    const t = nodeMap.get(d.target);
                    if (!s || !t) return "";
                    const mx = (s.x + t.x) / 2;
                    const my = (s.y + t.y) / 2;
                    return `M${s.x},${s.y} Q${mx},${my} ${t.x},${t.y}`;
                  });

                const radius = nodeRadius;

                const nodeSel = nodeGroup.selectAll("circle")
                  .data(nodes)
                  .enter()
                  .append("circle")
                  .attr("cx", (d) => d.x)
                  .attr("cy", (d) => d.y)
                  .attr("r", (d) => radius(d))
                  .attr("fill", (d) => d.color)
                  .attr("opacity", (d) => d.isSeed ? 1 : 0.76)
                  .attr("stroke", (d) => d.isSeed ? "#FFD166" : "none")
                  .attr("stroke-width", (d) => d.isSeed ? 2 : 0)
                  .style("cursor", "pointer");

                // Seed labels removed — tooltip on hover and side panel on click
                // provide node names without overlap issues at this density.
                // See Fruchterman & Reingold (1991) on label occlusion in dense layouts.

                const statsHtml = `
                  <div class="pt-overlay-title">Graph view</div>
                  <div>${nodes.length} nodes · ${edges.length} edges</div>
                  <div>${state.viewLabel || ""} · ${state.sortLabel || ""}</div>
                  <div style="margin-top:4px;color:#8aa4ff;">Largest component: ${state.globalMetrics?.largest_component_pct ?? "?"}%</div>
                `;
                statsBox.innerHTML = statsHtml;

                const kindMap = new Map();
                for (const _node of nodes) {
                  if (!kindMap.has(_node.kind)) kindMap.set(_node.kind, _node.color);
                }
                const edgeKindMap = new Map();
                for (const _edge of edges) {
                  if (!edgeKindMap.has(_edge.kind)) edgeKindMap.set(_edge.kind, _edge.color || "#bdbdbd");
                }
                legendBox.innerHTML = `<div class="pt-overlay-title">Legend</div>` +
                  `<div style="color:#7d8596;margin-bottom:4px;">Node types</div>` +
                  Array.from(kindMap.entries()).map(([k, c]) => `<div class="pt-legend-row"><span class="pt-swatch" style="background:${c}"></span>${k}</div>`).join("") +
                  `<div style="color:#7d8596;margin:8px 0 4px;">Edge types</div>` +
                  Array.from(edgeKindMap.entries()).map(([k, c]) => `<div class="pt-legend-row"><span class="pt-swatch" style="background:${c};border-radius:0;"></span>${k}</div>`).join("");

                leadsBox.innerHTML = leads.map((d) => `
                    <div class="pt-lead-card" data-id="${d.id}">
                      <div class="pt-lead-name" style="color:${d.color};">${d.name}</div>
                      <div class="pt-lead-detail">${d.kind} · score ${d.discoveryScore.toFixed(1)}</div>
                    </div>
                  `).join("");

                let selectedId = null;

                function resetView() {
                  selectedId = null;
                  edgeSel.attr("opacity", 0.10).attr("stroke-width", 0.6);
                  nodeSel.attr("opacity", (d) => d.isSeed ? 1 : 0.76).attr("r", (d) => radius(d));
                  sideBody.innerHTML = `<div class="pt-empty">Select a node to begin.</div>`;
                }

                function neighborInfo(nodeId) {
                  const connected = adj.get(nodeId) || [];
                  const neighbors = [];
                  const groups = new Map();
                  for (const _edge of connected) {
                    const otherId = _edge.source === nodeId ? _edge.target : _edge.source;
                    const other = nodeMap.get(otherId);
                    if (!other) continue;
                    neighbors.push(other);
                    if (!groups.has(_edge.kind)) groups.set(_edge.kind, []);
                    groups.get(_edge.kind).push(other);
                  }
                  return { connected, neighbors, groups };
                }

                function topBridgePaths(nodeId) {
                  const oneHop = adj.get(nodeId) || [];
                  const scored = [];
                  for (const e1 of oneHop) {
                    const midId = e1.source === nodeId ? e1.target : e1.source;
                    const midNode = nodeMap.get(midId);
                    if (!midNode) continue;
                    const secondHop = adj.get(midId) || [];
                    for (const e2 of secondHop) {
                      const dstId = e2.source === midId ? e2.target : e2.source;
                      if (dstId === nodeId) continue;
                      const dstNode = nodeMap.get(dstId);
                      if (!dstNode) continue;
                      const layerSpread = new Set([nodeMap.get(nodeId)?.layer, midNode.layer, dstNode.layer]).size;
                      const score = (midNode.crossLayerCount || 0) + (dstNode.discoveryScore || 0) + layerSpread * 4;
                      scored.push({ midNode, dstNode, e1, e2, score });
                    }
                  }
                  scored.sort((a, b) => b.score - a.score);
                  const seen = new Set();
                  const out = [];
                  for (const item of scored) {
                    const key = `${item.midNode.id}|${item.dstNode.id}|${item.e1.kind}|${item.e2.kind}`;
                    if (seen.has(key)) continue;
                    seen.add(key);
                    out.push(item);
                    if (out.length >= 6) break;
                  }
                  return out;
                }

                function highlight(nodeId) {
                  const current = nodeMap.get(nodeId);
                  if (!current) return;
                  selectedId = nodeId;
                  const connected = adj.get(nodeId) || [];
                  const neighborIds = new Set();
                  for (const _edge of connected) {
                    neighborIds.add(_edge.source === nodeId ? _edge.target : _edge.source);
                  }

                  edgeSel.attr("opacity", (d) => (d.source === nodeId || d.target === nodeId) ? 0.75 : 0.02)
                    .attr("stroke-width", (d) => (d.source === nodeId || d.target === nodeId) ? 1.8 : 0.4);
                  nodeSel.attr("opacity", (d) => {
                    if (d.id === nodeId) return 1;
                    return neighborIds.has(d.id) ? 0.95 : 0.06;
                  }).attr("r", (d) => {
                    if (d.id === nodeId) return 10;
                    if (neighborIds.has(d.id)) return Math.max(4, radius(d) + 1.2);
                    return radius(d);
                  });

                  const info = neighborInfo(nodeId);
                  const relationGroups = Array.from(info.groups.entries())
                    .sort((a, b) => b[1].length - a[1].length)
                    .map(([rel, arr]) => {
                      const chips = arr
                        .sort((a, b) => b.discoveryScore - a.discoveryScore)
                        .slice(0, 14)
                        .map((n) => `<span class="pt-chip" data-id="${n.id}">${n.name}</span>`)
                        .join("");
                      return `<div class="pt-card"><div style="font-weight:700;color:${state.edgeColors?.[rel] || '#8aa4ff'};margin-bottom:6px;">${rel} (${arr.length})</div><div class="pt-chip-wrap">${chips || '<span class="pt-empty">No nodes</span>'}</div></div>`;
                    })
                    .join("");

                  const paths = topBridgePaths(nodeId).map((p) => `
                    <div class="pt-card">
                      <div><span style="color:${p.midNode.color};font-weight:700;">${p.midNode.name}</span> → <span style="color:${p.dstNode.color};font-weight:700;">${p.dstNode.name}</span></div>
                      <div style="font-size:11px;color:#9aa5b5;margin-top:4px;">${current.kind} — ${p.e1.kind} → ${p.midNode.kind} — ${p.e2.kind} → ${p.dstNode.kind}</div>
                    </div>`).join("");

                  sideBody.innerHTML = `
                    <div class="pt-section">
                      <div class="pt-section-title">Overview</div>
                      <div class="pt-card">
                        <div style="font-size:18px;font-weight:700;color:${current.color};">${current.name}</div>
                        <div class="pt-sub">${current.kind} · ${current.layerName}</div>
                        <div class="pt-stat-grid">
                          <div class="pt-mini"><div class="v">${current.degree}</div><div class="l">degree</div></div>
                          <div class="pt-mini"><div class="v">${current.crossLayerCount}</div><div class="l">cross-layer</div></div>
                          <div class="pt-mini"><div class="v">${current.discoveryScore.toFixed(1)}</div><div class="l">score</div></div>
                        </div>
                      </div>
                    </div>
                    <div class="pt-section">
                      <div class="pt-section-title">Discover</div>
                      ${paths || '<div class="pt-empty">No 2-hop bridge paths found.</div>'}
                    </div>
                    <div class="pt-section">
                      <div class="pt-section-title">Neighbors by relation</div>
                      ${relationGroups || '<div class="pt-empty">No visible neighbors in current preset.</div>'}
                    </div>
                  `;

                  sideBody.querySelectorAll(".pt-chip").forEach((chip) => {
                    chip.addEventListener("click", () => highlight(chip.getAttribute("data-id")));
                  });
                }

                nodeSel
                  .on("mouseover", function(event, d) {
                    const info = neighborInfo(d.id);
                    tooltip.style.display = "block";
                    tooltip.innerHTML = `
                      <div style="font-weight:700;color:${d.color};">${d.name}</div>
                      <div style="font-size:10px;color:#9aa5b5;">${d.kind} · ${d.layerName}</div>
                      <div style="margin-top:4px;">degree ${d.degree} · score ${d.discoveryScore.toFixed(1)} · neighbors ${info.neighbors.length}</div>
                    `;
                    tooltip.style.left = `${event.clientX + 12}px`;
                    tooltip.style.top = `${event.clientY - 18}px`;
                  })
                  .on("mousemove", function(event) {
                    tooltip.style.left = `${event.clientX + 12}px`;
                    tooltip.style.top = `${event.clientY - 18}px`;
                  })
                  .on("mouseout", function() {
                    tooltip.style.display = "none";
                  })
                  .on("click", function(event, d) {
                    event.stopPropagation();
                    if (selectedId === d.id) {
                      resetView();
                    } else {
                      highlight(d.id);
                    }
                  });

                svg.on("click", () => resetView());

                leadsBox.querySelectorAll(".pt-lead-card").forEach((item) => {
                  item.addEventListener("click", () => highlight(item.getAttribute("data-id")));
                });

                if (searchTerm) {
                  const match = nodes.find((n) => n.name.toLowerCase().includes(searchTerm));
                  if (match) {
                    highlight(match.id);
                    const t = d3.zoomIdentity.translate(width / 2 - match.x * 1.6, height / 2 - match.y * 1.6).scale(1.6);
                    svg.transition().duration(500).call(d3.zoom().scaleExtent([0.35, 7]).on("zoom", (event) => {
                      zoomLayer.attr("transform", event.transform);
                    }).transform, t);
                  }
                }
              };

              draw();
              model.on("change:payload", draw);
            }

            export default { render };
            '''
        
    return (GraphWidget,)


@app.cell
def _(GraphWidget, anywidget_install_help, mo, payload_json):
    if GraphWidget is None:
        graph_widget = anywidget_install_help
    else:
        graph_widget = mo.ui.anywidget(GraphWidget(payload=payload_json))
    return (graph_widget,)


@app.cell
def _(edge_preset, graph_widget, mo, search_input, sort_mode, summary_md):
    controls = mo.hstack([search_input, edge_preset, sort_mode], justify="start", gap=1)
    mo.vstack([summary_md, controls, graph_widget])
    return


if __name__ == "__main__":
    app.run()