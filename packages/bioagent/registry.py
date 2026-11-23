"""Skill registry.

Central registry of every available skill and its entry script. The CLI and the
validation scripts both read from here, so the registration data lives in one
place instead of being scattered across several files.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

__all__ = ["SkillSpec", "REGISTRY", "get_skill", "list_skills", "skill_names"]


@dataclass(frozen=True)
class SkillSpec:
    """Registration record for a single skill."""

    name: str
    category: str
    summary: str
    entry: str
    keywords: tuple[str, ...] = ()

    @property
    def entry_path(self) -> Path:
        """Absolute path of the entry script."""
        root = Path(__file__).resolve().parents[2]
        return root / "skills" / self.name / self.entry

    @property
    def skill_md(self) -> Path:
        """Absolute path of SKILL.md."""
        root = Path(__file__).resolve().parents[2]
        return root / "skills" / self.name / "SKILL.md"


REGISTRY: dict[str, SkillSpec] = {
    spec.name: spec
    for spec in [
        SkillSpec(
            name="bulk-rnaseq",
            category="Transcriptomics",
            summary="RNA-seq quantification and differential expression: FastQC QC, Salmon quantification, DESeq2 enrichment",
            entry="scripts/run_pipeline.py",
            keywords=("rna-seq", "salmon", "deseq2", "diffexp", "transcriptome"),
        ),
        SkillSpec(
            name="single-cell-rnaseq",
            category="Transcriptomics",
            summary="Standard single-cell RNA-seq workflow: QC, normalization, PCA, clustering, cell annotation",
            entry="scripts/run_pipeline.py",
            keywords=("single-cell", "scrna", "scanpy", "anndata", "clustering"),
        ),
        SkillSpec(
            name="protein-structure",
            category="Structural Biology",
            summary="PDB structure parsing and secondary structure analysis: DSSP annotation, structure statistics, UniProt metadata",
            entry="scripts/analyze_structure.py",
            keywords=("pdb", "dssp", "structure", "secondary-structure", "uniprot"),
        ),
        SkillSpec(
            name="protein-interactions",
            category="Structural Biology",
            summary="Protein-ligand interaction analysis: PLIP parsing, residue extraction, interaction network visualization",
            entry="scripts/analyze_interactions.py",
            keywords=("plip", "interaction", "hydrogen-bond", "ligand", "docking"),
        ),
        SkillSpec(
            name="molecular-dynamics",
            category="Structural Biology",
            summary="Molecular dynamics simulation and trajectory analysis: OpenMM simulation, RMSD/RMSF/PCA/hydrogen bond analysis",
            entry="scripts/run_analysis.py",
            keywords=("md", "openmm", "trajectory", "rmsd", "simulation"),
        ),
        SkillSpec(
            name="pymol-structure-visualization",
            category="Structural Biology",
            summary="Inter-chain residue contact analysis: grid spatial partitioning + residue chemistry classification, output hotspot residues and PyMOL scripts",
            entry="scripts/analyze_contacts.py",
            keywords=("pymol", "contact", "interface", "hydrogen-bond",
                     "salt-bridge", "hotspot", "visualization"),
        ),
        SkillSpec(
            name="peptide-design",
            category="Protein Engineering",
            # Note: this describes the **actually implemented** capability — a
            # self-built lightweight CNN classifier — not "PepMLM generation".
            # The old repo's design_pepmlm/ depends on a 404MB pretrained model
            # file and requires a GPU; this skill does not bundle that model,
            # so it cannot claim it.
            summary="Antimicrobial peptide analysis and screening: physicochemical properties, CNN activity prediction, sequence clustering for redundancy removal",
            entry="scripts/analyze_peptides.py",
            keywords=("peptide", "amp", "antimicrobial", "cnn", "clustering"),
        ),
        SkillSpec(
            name="literature-search",
            category="Literature Mining",
            summary="Academic literature search and screening: multi-source retrieval and deduplication across PubMed / OpenAlex / Semantic Scholar",
            entry="scripts/search_literature.py",
            keywords=("pubmed", "openalex", "literature", "citation", "papers"),
        ),
        SkillSpec(
            name="scientific-figures",
            category="Scientific Visualization",
            summary="Publication-grade figure generation: ROC curves, confusion matrices, heatmaps, volcano plots, SHAP summaries",
            entry="scripts/make_figures.py",
            keywords=("figure", "plot", "roc", "heatmap", "publication"),
        ),
        SkillSpec(
            name="ena-download",
            category="Transcriptomics",
            summary="ENA raw sequencing data download: Portal API FASTQ links, resumable transfers, MD5 verification",
            entry="scripts/download_runs.py",
            keywords=("ena", "sra", "fastq", "sequencing", "rnaseq"),
        ),
        SkillSpec(
            name="alphafold-db",
            category="Structural Biology",
            summary="AlphaFold DB structure download: API link resolution, batch retrieval of predicted structures and pLDDT confidence",
            entry="scripts/download_structures.py",
            keywords=("alphafold", "afdb", "structure", "plddt", "prediction"),
        ),
        SkillSpec(
            name="bioproject-search",
            category="Literature Mining",
            summary="NCBI BioProject metadata search: E-utilities title/description/species retrieval, keyword filtering",
            entry="scripts/search_bioprojects.py",
            keywords=("bioproject", "ncbi", "e-utilities", "prjna", "metadata"),
        ),
        SkillSpec(
            name="ml-peptide-generation",
            category="Protein Engineering",
            summary="ML peptide generation: PepMLM masked language model candidate sampling, pseudo-perplexity ranking, reproducible",
            entry="scripts/generate_peptides.py",
            keywords=("peptide", "pepmlm", "binder", "masked-lm", "generation"),
        ),
        SkillSpec(
            name="bio-database-fetch",
            category="Cheminformatics",
            summary="PubChem compound and NCBI Taxonomy batch lookup: properties, CAS numbers, species taxonomy",
            entry="scripts/fetch_compounds.py",
            keywords=("pubchem", "taxonomy", "compound", "cas", "smiles", "ncbi"),
        ),
        SkillSpec(
            name="protein-docking-hdock",
            category="Structural Biology",
            summary="HDOCK protein docking: HTTP job submission, polling for Docking Score and confidence, model download",
            entry="scripts/submit_docking.py",
            keywords=("docking", "hdock", "protein-protein", "complex", "score"),
        ),
        SkillSpec(
            name="protein-docking-lightdock",
            category="Structural Biology",
            summary="LightDock local protein docking: driving GSO swarm optimization, parsing GSO scores, building the best complex",
            entry="scripts/run_docking.py",
            keywords=("docking", "lightdock", "gso", "swarm", "dfire", "local"),
        ),
        SkillSpec(
            name="paper-fulltext-extract",
            category="Literature Mining",
            summary="Structured full-text extraction from PDF: GROBID to TEI, section tree/figure/reference parsing, corpus search, lead sentences",
            entry="scripts/extract_fulltext.py",
            keywords=(
                "pdf", "tei", "grobid", "fulltext", "pdf-parsing",
                "corpus", "literature-mining", "introduction",
            ),
        ),
    ]
}


def get_skill(name: str) -> Optional[SkillSpec]:
    """Look up a skill by name."""
    return REGISTRY.get(name)


def list_skills() -> list[SkillSpec]:
    """Return every registered skill."""
    return list(REGISTRY.values())


def skill_names() -> list[str]:
    """Return every registered skill name."""
    return list(REGISTRY.keys())
