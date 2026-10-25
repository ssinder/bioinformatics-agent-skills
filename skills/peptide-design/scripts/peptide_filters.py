#!/usr/bin/env python3
"""多肽分类判别规则库：基于 GO 术语、关键词与理化性质的多维打分。

提供三个层面的判别能力：

1. **注释层面** —— GO 术语与 UniProt keyword 的加权打分
   （直接抗菌 / 毒素 / 免疫 / 细胞毒 / 抗病毒 / 其他防御）
2. **理化层面** —— 基于序列性质的抗菌倾向打分
   （疏水性、电荷、长度、两性、芳香性）
3. **综合层面** —— 两层加权合成，用于"是否为抗菌肽"与
   "是否为非抗菌肽"的二分类判定

打分规则源自 UniProt 短肽数据筛选流程（``filter_non_antimicrobial_strict.py``
与 ``export_to_csv.py``）中沉淀下来的 GO/关键词词表。
"""

from __future__ import annotations

import re
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))

from peptide_properties import clean_sequence  # noqa: E402

# ---------------------------------------------------------------- GO 词表

# 直接抗菌活性相关 GO
DIRECT_ANTIMICROBIAL_GO = {
    "GO:0042742",  # defense response to bacterium
    "GO:0050830",  # defense response to Gram-positive bacterium
    "GO:0050832",  # defense response to fungus
    "GO:0050829",  # defense response to Gram-negative bacterium
    "GO:0045087",  # innate immune response
    "GO:0042745",  # defense response to protozoan
}

# 杀伤其他生物细胞
KILLING_GO = {
    "GO:0031640",  # killing of cells of another organism
    "GO:0051715",  # cellular killing
    "GO:0019875",  # cytolysis
}

# 毒素活性
TOXIN_GO = {
    "GO:0090729",  # toxin activity
    "GO:0044398",  # ion channel inhibitor activity
    "GO:0044487",  # toxin activity (legacy)
    "GO:0044478",  # antibacterial toxin activity
    "GO:0035738",  # pore-forming toxin activity
    "GO:0044562",  # channel activity (pore forming)
    "GO:0051673",  # hemolysis
    "GO:0001897",  # cytolysis by hemoglobin
}

# 免疫相关
IMMUNE_GO = {
    "GO:0006952",  # defense response
    "GO:0006955",  # immune response
    "GO:0044419",  # activation of immune response
}

# 抗病毒
VIRAL_GO = {
    "GO:0019048",  # viral process
    "GO:0019079",  # viral replication
}

# 细胞外定位
EXTRACELLULAR_GO = {
    "GO:0005576",  # extracellular region
    "GO:0005615",  # extracellular space
    "GO:0062023",  # collagen-containing extracellular matrix
}

# ---------------------------------------------------------------- 关键词词表

DIRECT_ANTIMICROBIAL_KEYWORDS = {
    "antimicrobial", "defensin", "antibiotic", "bacteriocin",
    "antifungal", "antiviral", "antiparasitic", "antimicrobial peptide",
    "antimicrobial peptides", "antibacterial", "antibacterial peptide",
}

TOXIN_KEYWORDS = {
    "toxin", "neurotoxin", "cardiotoxin", "hemorrhagic toxin",
    "myotoxin", "dermonecrotic toxin", "enterotoxin", "hemolysin",
    "cytotoxic", "cytotoxin", "venom", "poison",
}

IMMUNE_KEYWORDS = {
    "immunity", "innate immunity", "adaptive immunity", "defense",
    "antibiotic resistance", "immune response", "immunomodulatory",
}

OTHER_DEFENSE_KEYWORDS = {
    "plant defense", "amphibian defense peptide", "insect defense",
    "fungicide", "pesticide", "herbicicide", "wound healing",
}

VIRAL_KEYWORDS = {
    "virulence", "viral", "virus", "viral protein", "influenza",
    "human immunodeficiency virus", "hepatitis", "coronavirus",
}

# ---------------------------------------------------------------- 打分权重

GO_WEIGHTS = {
    "direct_antimicrobial": 5.0,
    "killing": 5.0,
    "toxin": 4.0,
    "immune": 3.0,
    "viral": 2.0,
}

KEYWORD_WEIGHTS = {
    "direct_antimicrobial": 4.0,
    "toxin": 3.0,
    "immune": 2.0,
    "other_defense": 2.0,
    "viral": 2.0,
}


def _contains_any(text: str, needles: set[str]) -> bool:
    """判断文本是否包含词表中任一词（大小写不敏感，子串匹配）。"""
    if not text:
        return False
    lowered = text.lower()
    return any(n.lower() in lowered for n in needles)


def extract_go_ids(text: str) -> set[str]:
    """从自由文本中抽取 GO 编号（``GO:0042742`` 形式）。"""
    if not text:
        return set()
    return set(re.findall(r"GO:\d{7}", str(text)))


# UniProt 的 GoTerm 属性以 ``F:`` / ``P:`` / ``C:`` 标注三大类别
_GO_ASPECT_PREFIX = {"F": "molecular_function", "P": "biological_process", "C": "cellular_component"}

# 常见 GO 编号的 aspect 速查表（UniProt 官方 aspect）。
# 当输入只有 ``GO:xxxxxxx: term name`` 而缺少 aspect 前缀时用于兜底还原。
GO_ASPECT_HINTS: dict[str, str] = {
    # 分子功能 (F)
    "GO:0003674": "molecular_function",
    "GO:0005488": "molecular_function",
    "GO:0003824": "molecular_function",
    "GO:0005515": "molecular_function",
    "GO:0003735": "molecular_function",
    "GO:0004252": "molecular_function",
    "GO:0004251": "molecular_function",
    "GO:0016301": "molecular_function",
    "GO:0090729": "molecular_function",
    "GO:0044487": "molecular_function",
    "GO:0044478": "molecular_function",
    "GO:0035738": "molecular_function",
    "GO:0044398": "molecular_function",
    "GO:0044562": "molecular_function",
    # 生物过程 (P)
    "GO:0006952": "biological_process",
    "GO:0006955": "biological_process",
    "GO:0045087": "biological_process",
    "GO:0044419": "biological_process",
    "GO:0042742": "biological_process",
    "GO:0042745": "biological_process",
    "GO:0050830": "biological_process",
    "GO:0050832": "biological_process",
    "GO:0050829": "biological_process",
    "GO:0031640": "biological_process",
    "GO:0051715": "biological_process",
    "GO:0019875": "biological_process",
    "GO:0001897": "biological_process",
    "GO:0051673": "biological_process",
    "GO:0019048": "biological_process",
    "GO:0019079": "biological_process",
    # 细胞组分 (C)
    "GO:0005576": "cellular_component",
    "GO:0005615": "cellular_component",
    "GO:0062023": "cellular_component",
}


def split_go_terms(text: str) -> dict[str, list[str]]:
    """把合并的 GO 文本按三大类别拆分。

    UniProt 返回的 GoTerm 属性形如 ``GO:0042742: F:defense response to bacterium``，
    其中 ``F:`` / ``P:`` / ``C:`` 分别代表分子功能 / 生物过程 / 细胞组分。
    若缺少 aspect 前缀，则回退到 ``GO_ASPECT_HINTS`` 速查表；
    两者都无法判定时归入 ``unknown``。

    Args:
        text: 以 ``;`` 分隔的 GO 术语文本。

    Returns:
        ``{"molecular_function": [...], "biological_process": [...],
        "cellular_component": [...], "unknown": [...]}``
    """
    buckets: dict[str, list[str]] = {name: [] for name in _GO_ASPECT_PREFIX.values()}
    buckets["unknown"] = []

    for chunk in re.split(r"[;,|]", str(text or "")):
        term = chunk.strip()
        if not term:
            continue
        aspect = None
        for prefix, name in _GO_ASPECT_PREFIX.items():
            if re.search(rf"\b{prefix}:", term):
                aspect = name
                break
        if aspect is None:
            go_ids = extract_go_ids(term)
            for go_id in go_ids:
                if go_id in GO_ASPECT_HINTS:
                    aspect = GO_ASPECT_HINTS[go_id]
                    break
        buckets[aspect or "unknown"].append(term)

    return buckets


def annotation_score(row: dict[str, Any]) -> dict[str, Any]:
    """基于 GO 术语与关键词计算注释层面得分。

    Args:
        row: 记录字典，需含 ``go_terms`` 与 ``keywords`` 字段（可缺省）。

    Returns:
        ``{"score": float, "categories": {...}, "hit_go": [...], "hit_keywords": [...]}``
    """
    go_text = str(row.get("go_terms", "") or "")
    kw_text = str(row.get("keywords", "") or "")
    go_ids = extract_go_ids(go_text)

    categories = {
        "direct_antimicrobial": False,
        "killing": False,
        "toxin": False,
        "immune": False,
        "viral": False,
        "extracellular": False,
    }
    score = 0.0

    if go_ids & DIRECT_ANTIMICROBIAL_GO:
        categories["direct_antimicrobial"] = True
        score += GO_WEIGHTS["direct_antimicrobial"]
    if go_ids & KILLING_GO:
        categories["killing"] = True
        score += GO_WEIGHTS["killing"]
    if go_ids & TOXIN_GO:
        categories["toxin"] = True
        score += GO_WEIGHTS["toxin"]
    if go_ids & IMMUNE_GO:
        categories["immune"] = True
        score += GO_WEIGHTS["immune"]
    if go_ids & VIRAL_GO:
        categories["viral"] = True
        score += GO_WEIGHTS["viral"]
    if go_ids & EXTRACELLULAR_GO:
        categories["extracellular"] = True

    kw_categories = {
        "direct_antimicrobial": _contains_any(kw_text, DIRECT_ANTIMICROBIAL_KEYWORDS),
        "toxin": _contains_any(kw_text, TOXIN_KEYWORDS),
        "immune": _contains_any(kw_text, IMMUNE_KEYWORDS),
        "other_defense": _contains_any(kw_text, OTHER_DEFENSE_KEYWORDS),
        "viral": _contains_any(kw_text, VIRAL_KEYWORDS),
    }
    for name, hit in kw_categories.items():
        if hit:
            score += KEYWORD_WEIGHTS[name]

    hit_keywords = [
        kw.strip()
        for kw in re.split(r"[;,|]", kw_text)
        if kw.strip()
        and (
            _contains_any(kw, DIRECT_ANTIMICROBIAL_KEYWORDS)
            or _contains_any(kw, TOXIN_KEYWORDS)
            or _contains_any(kw, IMMUNE_KEYWORDS)
            or _contains_any(kw, OTHER_DEFENSE_KEYWORDS)
            or _contains_any(kw, VIRAL_KEYWORDS)
        )
    ]

    return {
        "score": round(score, 2),
        "categories": categories,
        "keyword_categories": kw_categories,
        "hit_go": sorted(go_ids),
        "hit_keywords": hit_keywords,
    }


def physicochemical_score(props: dict[str, Any]) -> dict[str, Any]:
    """基于理化性质计算抗菌倾向得分。

    规则（各条独立计分后相加，范围约 0–10）：

    - 疏水性适中偏高（Kyte-Doolittle 1.0–4.0）：+2.0
    - 存在净正电荷（pH 7 下 > 0）：+2.0（膜结合的关键因素）
    - 长度在 12–50 aa 区间：+2.0（典型 AMP 长度）
    - 芳香性 0.05–0.25：+1.0
    - Boman 值 > 0：+1.5（潜在抗菌活性指数）
    - 含半胱氨酸（可能形成二硫键，多见于防御肽）：+0.5
    - 不稳定指数 < 40：+1.0
    """
    score = 0.0
    reasons: list[str] = []

    kd = float(props.get("kyte_doolittle", 0.0))
    if 1.0 <= kd <= 4.0:
        score += 2.0
        reasons.append(f"疏水性适中 (KD={kd:.2f})")
    elif kd > 4.0:
        score += 0.5
        reasons.append(f"疏水性强但可能溶解性差 (KD={kd:.2f})")

    charge = float(props.get("net_charge_ph7", 0.0))
    if charge > 0:
        score += 2.0
        reasons.append(f"净正电荷 {charge:.2f}")
    elif charge < -3:
        reasons.append(f"净负电荷 {charge:.2f}，不利于膜结合")

    length = int(props.get("length", 0))
    if 12 <= length <= 50:
        score += 2.0
        reasons.append(f"长度 {length} aa 落在典型 AMP 区间")
    elif length < 8:
        reasons.append(f"长度 {length} aa 过短")

    aromatic = float(props.get("aromaticity", 0.0))
    if 0.05 <= aromatic <= 0.25:
        score += 1.0
        reasons.append(f"芳香性适中 ({aromatic:.3f})")

    boman = float(props.get("boman_index", 0.0))
    if boman > 0:
        score += 1.5
        reasons.append(f"Boman 值 {boman:.2f} > 0")

    if int(props.get("cysteine_count", 0)) > 0:
        score += 0.5

    if float(props.get("instability_index", 100.0)) < 40:
        score += 1.0
        reasons.append("不稳定指数 < 40")

    return {"score": round(min(score, 10.0), 2), "reasons": reasons}


def combined_score(
    annotation: dict[str, Any],
    physchem: dict[str, Any],
    annotation_weight: float = 0.6,
) -> dict[str, Any]:
    """合成注释层与理化层得分。

    ``total = annotation_weight × annotation_score + (1-w) × physchem_score``

    归一化到 0–10 量纲，便于与阈值比较。
    """
    total = (
        annotation_weight * float(annotation.get("score", 0.0))
        + (1.0 - annotation_weight) * float(physchem.get("score", 0.0))
    )
    return {
        "total": round(total, 3),
        "annotation_score": annotation.get("score", 0.0),
        "physicochemical_score": physchem.get("score", 0.0),
    }


def classify(
    row: dict[str, Any],
    props: dict[str, Any],
    antimicrobial_threshold: float = 5.0,
    non_antimicrobial_threshold: float = 2.0,
    annotation_weight: float = 0.6,
) -> dict[str, Any]:
    """对一条肽做三分类判定。

    Returns:
        ``{"label": "antimicrobial"|"non_antimicrobial"|"uncertain", "total": float, ...}``
    """
    annotation = annotation_score(row)
    physchem = physicochemical_score(props)
    combined = combined_score(annotation, physchem, annotation_weight)

    total = combined["total"]
    if total >= antimicrobial_threshold:
        label = "antimicrobial"
    elif total <= non_antimicrobial_threshold:
        label = "non_antimicrobial"
    else:
        label = "uncertain"

    return {
        "label": label,
        "total": total,
        "annotation_score": combined["annotation_score"],
        "physicochemical_score": combined["physicochemical_score"],
        "categories": annotation["categories"],
        "hit_keywords": annotation["hit_keywords"],
        "reasons": physchem["reasons"],
    }


def quality_score(row: dict[str, Any], props: dict[str, Any]) -> float:
    """样本质量综合分（用于"高质量样本提取"）。

    组成：

    - UniProt 蛋白质存在水平（Evidence at protein level = 最高档）
    - 序列完整（无未知残基 X）
    - 注释完整度（GO 术语与关键词的条目数）
    - 理化性质在合理区间（不稳定指数 < 40、长度合理）
    """
    score = 0.0

    existence = str(row.get("protein_existence", "") or "")
    if "Evidence at protein level" in existence:
        score += 40
    elif "Evidence at transcript level" in existence:
        score += 25
    elif "Protein inferred from homology" in existence:
        score += 15
    elif existence:
        score += 5

    if not props.get("has_unknown"):
        score += 20
    else:
        score -= 10

    go_count = len(extract_go_ids(str(row.get("go_terms", "") or "")))
    kw_count = len([k for k in re.split(r"[;,|]", str(row.get("keywords", "") or "")) if k.strip()])
    score += min(20, go_count * 2 + kw_count)

    length = int(props.get("length", 0))
    if 8 <= length <= 60:
        score += 10
    if float(props.get("instability_index", 100.0)) < 40:
        score += 10

    return round(score, 2)


def sequence_complexity(sequence: str) -> float:
    """序列复杂度（Shannon 熵，归一化到 0–1）。

    值越高说明序列组成越多样（Shannon entropy / log2(20)）。
    """
    sequence = clean_sequence(sequence)
    import math

    if not sequence:
        return 0.0
    counts: dict[str, int] = {}
    for ch in sequence:
        counts[ch] = counts.get(ch, 0) + 1
    total = len(sequence)
    entropy = 0.0
    for count in counts.values():
        p = count / total
        entropy -= p * math.log2(p)
    return round(entropy / math.log2(20), 4)


def gc_content(sequence: str) -> float:
    """GC 含量。"""
    sequence = clean_sequence(sequence)
    if not sequence:
        return 0.0
    gc = sum(1 for ch in sequence if ch in "GC")
    return round(gc / len(sequence), 4)


def cys_rich(sequence: str, threshold: float = 0.15) -> bool:
    """是否富含半胱氨酸（阈值按序列长度的比例）。"""
    sequence = clean_sequence(sequence)
    if not sequence:
        return False
    return sequence.count("C") / len(sequence) > threshold


def find_repeated_motifs(sequence: str, motif_length: int = 3, min_repeat: int = 3) -> list[str]:
    """查找重复出现超过阈值的短 motif（用于识别保守结构域）。"""
    sequence = clean_sequence(sequence)
    if len(sequence) < motif_length:
        return []
    counts: dict[str, int] = {}
    for i in range(len(sequence) - motif_length + 1):
        motif = sequence[i : i + motif_length]
        counts[motif] = counts.get(motif, 0) + 1
    return sorted(m for m, c in counts.items() if c >= min_repeat)


def annotate_go_distribution(records: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """统计 GO 三大类别的覆盖度（分子功能 / 生物过程 / 细胞组分）。

    优先读取已拆分的三列；若只有合并的 ``go_terms``，则按 UniProt 的
    ``F:`` / ``P:`` / ``C:`` aspect 前缀还原类别。
    """
    mf = bp = cc = 0
    total = len(records)
    for row in records:
        has_mf = str(row.get("go_molecular_function", "") or "").strip()
        has_bp = str(row.get("go_biological_process", "") or "").strip()
        has_cc = str(row.get("go_cellular_component", "") or "").strip()
        if not (has_mf or has_bp or has_cc):
            buckets = split_go_terms(str(row.get("go_terms", "") or ""))
            has_mf = bool(buckets["molecular_function"])
            has_bp = bool(buckets["biological_process"])
            has_cc = bool(buckets["cellular_component"])
        mf += 1 if has_mf else 0
        bp += 1 if has_bp else 0
        cc += 1 if has_cc else 0
    return {
        "total": total,
        "molecular_function": mf,
        "biological_process": bp,
        "cellular_component": cc,
        "molecular_function_pct": round(mf / total * 100, 2) if total else 0.0,
        "biological_process_pct": round(bp / total * 100, 2) if total else 0.0,
        "cellular_component_pct": round(cc / total * 100, 2) if total else 0.0,
    }
