"""技能注册表。

集中登记所有可用技能及其入口脚本，CLI 与校验脚本都从这里取数据，
避免注册信息分散在多处导致遗漏。
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

__all__ = ["SkillSpec", "REGISTRY", "get_skill", "list_skills", "skill_names"]


@dataclass(frozen=True)
class SkillSpec:
    """一个技能的注册信息。"""

    name: str
    category: str
    summary: str
    entry: str
    keywords: tuple[str, ...] = ()

    @property
    def entry_path(self) -> Path:
        """入口脚本的绝对路径。"""
        root = Path(__file__).resolve().parents[2]
        return root / "skills" / self.name / self.entry

    @property
    def skill_md(self) -> Path:
        """SKILL.md 的绝对路径。"""
        root = Path(__file__).resolve().parents[2]
        return root / "skills" / self.name / "SKILL.md"


REGISTRY: dict[str, SkillSpec] = {
    spec.name: spec
    for spec in [
        SkillSpec(
            name="bulk-rnaseq",
            category="转录组学",
            summary="RNA-seq 定量与差异表达分析：FastQC 质控、Salmon 定量、DESeq2 富集",
            entry="scripts/run_pipeline.py",
            keywords=("rna-seq", "salmon", "deseq2", "diffexp", "transcriptome"),
        ),
        SkillSpec(
            name="single-cell-rnaseq",
            category="转录组学",
            summary="单细胞 RNA-seq 标准流程：QC、归一化、PCA、聚类、细胞注释",
            entry="scripts/run_pipeline.py",
            keywords=("single-cell", "scrna", "scanpy", "anndata", "clustering"),
        ),
        SkillSpec(
            name="protein-structure",
            category="结构生物学",
            summary="PDB 结构解析与二级结构分析：DSSP 注释、结构统计、UniProt 元数据",
            entry="scripts/analyze_structure.py",
            keywords=("pdb", "dssp", "structure", "secondary-structure", "uniprot"),
        ),
        SkillSpec(
            name="protein-interactions",
            category="结构生物学",
            summary="蛋白质-配体相互作用分析：PLIP 解析、残基抽取、交互网络可视化",
            entry="scripts/analyze_interactions.py",
            keywords=("plip", "interaction", "hydrogen-bond", "ligand", "docking"),
        ),
        SkillSpec(
            name="molecular-dynamics",
            category="结构生物学",
            summary="分子动力学模拟与轨迹分析：OpenMM 模拟、RMSD/RMSF/PCA/氢键分析",
            entry="scripts/run_analysis.py",
            keywords=("md", "openmm", "trajectory", "rmsd", "simulation"),
        ),
        SkillSpec(
            name="pymol-structure-visualization",
            category="结构生物学",
            summary="链间残基接触分析：网格空间划分 + 残基化学分类，输出热点残基与 PyMOL 脚本",
            entry="scripts/analyze_contacts.py",
            keywords=("pymol", "contact", "interface", "hydrogen-bond",
                     "salt-bridge", "hotspot", "visualization"),
        ),
        SkillSpec(
            name="peptide-design",
            category="蛋白工程",
            # 注意：这里描述的是**实际实现的能力**——自研轻量 CNN 分类器，
            # 而不是「PepMLM 生成」。老仓库的 design_pepmlm/ 依赖 404MB 的
            # 预训练模型文件且需要 GPU，本技能并未内置该模型，因此不能宣称。
            summary="抗菌肽分析与筛选：理化性质计算、CNN 活性判别、序列聚类去冗余",
            entry="scripts/analyze_peptides.py",
            keywords=("peptide", "amp", "antimicrobial", "cnn", "clustering"),
        ),
        SkillSpec(
            name="literature-search",
            category="文献挖掘",
            summary="学术文献检索与筛选：PubMed / OpenAlex / Semantic Scholar 多源检索与去重",
            entry="scripts/search_literature.py",
            keywords=("pubmed", "openalex", "literature", "citation", "papers"),
        ),
        SkillSpec(
            name="scientific-figures",
            category="数据可视化",
            summary="顶刊级图表生成：ROC 曲线、混淆矩阵、热图、火山图、SHAP 摘要",
            entry="scripts/make_figures.py",
            keywords=("figure", "plot", "roc", "heatmap", "publication"),
        ),
        SkillSpec(
            name="ena-download",
            category="转录组学",
            summary="ENA 原始测序数据下载：Portal API 取 FASTQ 链接、断点续传、MD5 校验",
            entry="scripts/download_runs.py",
            keywords=("ena", "sra", "fastq", "sequencing", "rnaseq"),
        ),
        SkillSpec(
            name="alphafold-db",
            category="结构生物学",
            summary="AlphaFold DB 结构下载：API 解析下载链接、批量获取预测结构与 pLDDT 置信度",
            entry="scripts/download_structures.py",
            keywords=("alphafold", "afdb", "structure", "plddt", "prediction"),
        ),
        SkillSpec(
            name="bioproject-search",
            category="文献挖掘",
            summary="NCBI BioProject 元数据检索：E-utilities 取标题/描述/物种，支持关键词筛选",
            entry="scripts/search_bioprojects.py",
            keywords=("bioproject", "ncbi", "e-utilities", "prjna", "metadata"),
        ),
        SkillSpec(
            name="ml-peptide-generation",
            category="蛋白工程",
            summary="ML 肽链生成：PepMLM 掩码语言模型采样候选肽，伪困惑度排序，可复现",
            entry="scripts/generate_peptides.py",
            keywords=("peptide", "pepmlm", "binder", "masked-lm", "generation"),
        ),
        SkillSpec(
            name="bio-database-fetch",
            category="化学信息学",
            summary="PubChem 化合物与 NCBI Taxonomy 批量检索：属性、CAS 号、物种分类",
            entry="scripts/fetch_compounds.py",
            keywords=("pubchem", "taxonomy", "compound", "cas", "smiles", "ncbi"),
        ),
        SkillSpec(
            name="protein-docking-hdock",
            category="结构生物学",
            summary="HDOCK 蛋白质对接：HTTP 提交作业、轮询取 Docking Score 与置信度、下载模型",
            entry="scripts/submit_docking.py",
            keywords=("docking", "hdock", "protein-protein", "complex", "score"),
        ),
        SkillSpec(
            name="protein-docking-lightdock",
            category="结构生物学",
            summary="LightDock 本地蛋白质对接：驱动 GSO 群体优化、解析 GSO 打分、生成最优复合物",
            entry="scripts/run_docking.py",
            keywords=("docking", "lightdock", "gso", "swarm", "dfire", "local"),
        ),
        SkillSpec(
            name="paper-fulltext-extract",
            category="文献挖掘",
            summary="PDF 文献全文结构化抽取：GROBID 转 TEI、章节树/图表/参考文献解析、语料检索、引言首句",
            entry="scripts/extract_fulltext.py",
            keywords=(
                "pdf", "tei", "grobid", "fulltext", "pdf-parsing",
                "corpus", "literature-mining", "introduction",
            ),
        ),
    ]
}


def get_skill(name: str) -> Optional[SkillSpec]:
    """按名称查找技能。"""
    return REGISTRY.get(name)


def list_skills() -> list[SkillSpec]:
    """返回全部已注册技能。"""
    return list(REGISTRY.values())


def skill_names() -> list[str]:
    """返回全部技能名。"""
    return list(REGISTRY.keys())
