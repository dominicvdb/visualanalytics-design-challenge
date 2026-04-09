# Biomedical Knowledge Graph Explorer

An interactive visual analytics interface for exploring a subset of [Hetionet](https://het.io/), a biomedical knowledge graph with 47,031 nodes and 2.25 million relationships assembled from 29 public databases. Built with [Marimo](https://marimo.io/) and [D3.js](https://d3js.org/).

The interface uses **semantic substrates** (Shneiderman & Aris, 2006) to organize 11 entity types into four horizontal layers — Clinical, Pharmaceutical, Genomic, and Biological — with force-directed node placement within each layer (Fruchterman & Reingold, 1991).

## Features

- **Layered semantic substrates**: Nodes are grouped into four biomedical layers with force-directed overlap removal within each layer
- **Focus+context interaction**: Click any node to highlight its connections while the rest of the graph dims as context
- **Discovery panel**: Side panel showing node overview, two-hop bridge paths, and neighbors grouped by relation type
- **Discovery leads**: Top-scored nodes shown as suggested starting points for exploration
- **Search**: Type a node name to zoom directly to it
- **Relationship presets**: Filter visible edge types by analytical perspective (cross-layer discovery, clinical focus, biological mechanisms)
- **Sorting**: Reorder nodes within layers by degree, discovery score, alphabetical, or seed proximity

## Setup

### Requirements

- Python 3.10+
- A modern browser (Chrome, Firefox, Edge)

### Installation

```bash
git clone https://github.com/YOUR_USERNAME/YOUR_REPO.git
cd YOUR_REPO
pip install -r requirements.txt
```

### Running

```bash
marimo run pathtrace_app.py
```

The app will open in your browser.

## Project structure

```
├── pathtrace_app.py          # Marimo app (main file)
├── pyproject.toml             # Marimo runtime config
├── requirements.txt           # Python dependencies
├── subset_hetionet.py         # Script to regenerate the subset from full Hetionet
├── data/
│   └── hetionet_subset.json   # Pre-built subset (~1,700 nodes, ~5,000 edges)
└── README.md
```

## Data

The app uses a disease-focused subset of Hetionet v1.0. The subset was built by:

1. Selecting 11 seed diseases (Alzheimer's, breast cancer, diabetes, etc.) that have known drug treatments
2. Performing a 2-hop BFS from those seeds to pull in connected compounds, genes, pathways, and other entities
3. Capping each entity type to keep the graph under 3,000 nodes
4. Applying edge budgets per relationship type for browser performance

To regenerate the subset from the full Hetionet dataset:

1. Download `hetionet-v1.0.json.bz2` from the [Hetionet GitHub repository](https://github.com/hetio/hetionet)
2. Place it in `data/`
3. Run `python subset_hetionet.py`

## Dataset citation

Himmelstein, D. S., Lizee, A., Hessler, C., Brueggeman, L., Chen, S. L., Hadley, D., Green, A., Khankhanian, P., & Baranzini, S. E. (2017). Systematic integration of biomedical knowledge prioritizes drugs for repurposing. *eLife*, *6*, e26726. https://doi.org/10.7554/eLife.26726
