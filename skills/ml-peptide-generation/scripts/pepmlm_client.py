#!/usr/bin/env python3
"""PepMLM 掩码语言模型客户端。

对应老仓库 ``design_pepmlm/scripts/`` 下的三个生成脚本
（``batch_generate_peptides.py`` / ``test_pepmlm.py`` /
``test_new_structure.py``）。这三个脚本共享同一套核心逻辑，本模块把它
抽出来，并修掉一批真实缺陷。

老脚本的共同缺陷
----------------
1. **模型在 import 时就加载**。``test_pepmlm.py`` 第 18-24 行在模块顶层
   调用 ``AutoModelForMaskedLM.from_pretrained``——只要 import 这个文件，
   就会去 HuggingFace 拉 650M 参数的模型。``--help`` 因此也无法执行，
   而"import 一个模块不该有副作用"是本仓库的硬性约束。
2. **没有任何随机种子控制**。生成用 ``Categorical(probabilities).sample()``
   采样，同一输入跑两次得到不同结果，实验无法复现。
3. **硬编码 macOS 绝对路径**（``/Users/jung/Desktop/trae/tool/...``）。
4. **三个脚本近乎逐行重复**——``test_new_structure.py`` 与
   ``batch_generate_peptides.py`` 只差输入输出路径。
5. **Top-K 采样会把非法氨基酸选进来**。模型词表包含 ``X``（未知）、
   ``B``（ambiguous）这类特殊 token，直接采样会产生合成上无意义的肽链；
   老脚本只是在事后用 ``Has_Unknown_AA`` 标记一下，没有在采样时排除。
6. **"显示前 3 条结果"的循环变量越界**。``batch_generate_peptides.py``
   第 252 行在长度循环**之外**引用 ``results``，拿到的是最后一次长度
   的结果——当某个长度被 ``continue`` 跳过时会直接 ``UnboundLocalError``。

本模块的改进
------------
- 模型**惰性加载**：首次真正需要时才从本地缓存或 HuggingFace 拉取。
- **完整的种子控制**：``set_seed()`` 同时设置 ``random``、``numpy``、
  ``torch``（CPU 与 CUDA），并在模块内提供可注入的采样函数，使同样的
  种子必然产出同样的肽链。
- **采样时排除非法 token**：从 Top-K 候选里剔除非标准氨基酸，再在
  剩余候选中按概率采样；候选被全部剔除时才回退到原始分布并记录警告。
- **纯 Python 部分不依赖 torch**：FASTA 解析、结果筛选、参数统计
  都不碰深度学习栈，因此没有 GPU 的机器也能完成分析环节。

PepMLM 的工作方式
------------------
PepMLM 是一个**掩码语言模型**：把目标蛋白序列后面接 N 个 ``<mask>``，
模型预测每个 mask 位置上的氨基酸分布。生成 N 个氨基酸的肽链，就是
迭代 N 次、每次采样一个位置。

伪困惑度（pseudo perplexity）
-----------------------------
生成完成后，把「蛋白 + 肽链」一起送进模型，逐个 mask 掉肽链上的位置，
计算模型恢复这些位置的交叉熵损失，再取指数：

.. code-block:: text

    pseudo_perplexity = exp(mean_cross_entropy)

它衡量的是「模型认为这条肽链有多像它自己会生成的」。**数值越低越好**，
可以作为候选肽链的排序依据。注意它衡量的是模型自身的信念，而不是实验
结合亲和力——低 PPL 只意味着「模型觉得这条肽链合理」。
"""

from __future__ import annotations

import math
import random
import sys
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))

from biopipe_core import get_logger  # noqa: E402

__all__ = [
    "PepMLMConfig",
    "PepMLMClient",
    "set_seed",
    "VALID_AMINO_ACIDS",
    "PepMLMUnavailableError",
    "DEFAULT_MODEL_NAME",
]

LOGGER = get_logger("biopipe.ml-peptide-generation")

#: 默认的预训练模型
DEFAULT_MODEL_NAME = "TianlaiChen/PepMLM-650M"

#: 20 种标准氨基酸。模型词表里还有 X（未知）、B/Z（ambiguous）等特殊
#: token，采样时必须排除——见 PepMLMClient._filter_candidates。
VALID_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")

#: PepMLM 的最大上下文长度（token 数）。超出会被静默截断，导致蛋白
#: 序列被切掉一截而模型毫无察觉。
MAX_CONTEXT_LENGTH = 552


class PepMLMUnavailableError(RuntimeError):
    """torch / transformers 缺失，或模型无法加载。"""


def set_seed(seed: int) -> None:
    """设置全局随机种子。

    必须同时设置三个来源：Python 的 ``random``（序列打乱等）、
    ``numpy``（如果有）和 ``torch``（模型采样）。只设其中一个，
    结果仍然不可复现。

    Args:
        seed: 种子值。相同种子必然产出相同结果。
    """
    random.seed(seed)

    try:
        import numpy as np
    except ImportError:
        pass
    else:
        np.random.seed(seed % (2**32))

    try:
        import torch
    except ImportError:
        # 没装 torch 时不影响纯 Python 部分的复现性
        return

    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


@dataclass
class PepMLMConfig:
    """PepMLM 运行配置。"""

    model_name: str = DEFAULT_MODEL_NAME
    peptide_lengths: tuple[int, ...] = (10, 12, 15)
    top_k: int = 3
    num_binders: int = 10
    seed: int = 42
    device: Optional[str] = None
    batch_size: int = 1
    max_context_length: int = MAX_CONTEXT_LENGTH
    exclude_invalid_amino_acids: bool = True
    extra: dict[str, Any] = field(default_factory=dict)


class PepMLMClient:
    """PepMLM 掩码语言模型的薄封装。

    模型**惰性加载**：只有第一次真正需要推理时才加载，因此构造这个
    对象本身是廉价的，``--help`` 与参数校验都能在没装 torch 的机器上跑。
    """

    def __init__(self, config: Optional[PepMLMConfig] = None) -> None:
        self.config = config or PepMLMConfig()
        self._tokenizer = None
        self._model = None
        self._device = None
        self._valid_token_ids: Optional[set[int]] = None
        self._rng = random.Random(self.config.seed)

    # ------------------------------------------------------------------
    # 模型加载
    # ------------------------------------------------------------------
    def _ensure_model(self) -> None:
        """惰性加载模型与 tokenizer。"""
        if self._model is not None:
            return

        try:
            import torch
            from transformers import AutoModelForMaskedLM, AutoTokenizer
        except ImportError as exc:
            raise PepMLMUnavailableError(
                "生成肽链需要 torch 与 transformers。\n"
                "请先安装: pip install torch transformers\n"
                "（若只想做结果筛选与参数分析，可直接用 analyze_peptides.py，"
                "它不依赖深度学习栈）"
            ) from exc

        device = self.config.device or ("cuda" if torch.cuda.is_available() else "cpu")
        LOGGER.info("正在加载模型 %s 到 %s …", self.config.model_name, device)

        tokenizer = AutoTokenizer.from_pretrained(self.config.model_name)
        model = AutoModelForMaskedLM.from_pretrained(self.config.model_name)
        model = model.to(device)
        model.eval()

        self._tokenizer = tokenizer
        self._model = model
        self._device = device
        self._build_valid_token_ids()

    def _build_valid_token_ids(self) -> None:
        """收集「标准氨基酸」对应的 token id 集合。

        模型的词表里混有 ``<mask>``、``<cls>``、``<pad>`` 这些特殊
        token，以及 ``X``、``B``、``Z`` 这些 ambigous/未知氨基酸。
        直接在完整词表上做 Top-K 采样，命中它们的概率不低，生成的肽链
        会含有 ``X`` 这类合成上无意义的残基。

        这里预先算出「合法 token id」集合，采样时用掩码把非法候选
        排除掉——比生成完再过滤更有效，因为它不会浪费采样次数。
        """
        tokenizer = self._tokenizer
        valid: set[int] = set()
        for amino_acid in sorted(VALID_AMINO_ACIDS):
            token_id = tokenizer.convert_tokens_to_ids(amino_acid)
            if token_id is not None and token_id != tokenizer.unk_token_id:
                valid.add(int(token_id))
        self._valid_token_ids = valid
        LOGGER.debug("合法氨基酸 token 数: %d", len(valid))

    # ------------------------------------------------------------------
    # 生成
    # ------------------------------------------------------------------
    def generate_peptides(
        self,
        protein_seq: str,
        peptide_length: int = 15,
        top_k: Optional[int] = None,
        num_binders: Optional[int] = None,
    ) -> list[dict]:
        """为一条蛋白序列生成候选肽链。

        Args:
            protein_seq: 靶标蛋白序列（标准氨基酸）。
            peptide_length: 肽链长度。
            top_k: 每个 mask 位置保留的候选数。
            num_binders: 生成多少条候选。

        Returns:
            每条候选一个字典，含 ``binder``、``pseudo_perplexity``、
            ``binder_length``、``has_unknown_aa``。
        """
        # 必须先确保模型可用——它负责检查 torch/transformers 是否安装，
        # 缺失时抛带安装指引的 PepMLMUnavailableError，而不是裸的
        # ModuleNotFoundError。
        self._ensure_model()
        import torch

        peptide_length = int(peptide_length)
        top_k = int(top_k if top_k is not None else self.config.top_k)
        num_binders = int(num_binders if num_binders is not None else self.config.num_binders)

        if peptide_length <= 0:
            raise ValueError(f"肽链长度必须为正整数，收到 {peptide_length}")
        if top_k <= 0:
            raise ValueError(f"top_k 必须为正整数，收到 {top_k}")
        if num_binders <= 0:
            raise ValueError(f"生成数量必须为正整数，收到 {num_binders}")

        total_length = len(protein_seq) + peptide_length
        if total_length > self.config.max_context_length:
            raise ValueError(
                f"蛋白长度 ({len(protein_seq)}) + 肽链长度 ({peptide_length}) = "
                f"{total_length} 超过模型上下文上限 {self.config.max_context_length}。"
                "请改用更短的靶标片段。"
            )

        results: list[dict] = []
        for _ in range(num_binders):
            binder = self._sample_one(protein_seq, peptide_length, top_k, torch)
            ppl = self.compute_pseudo_perplexity(protein_seq, binder)
            results.append(
                {
                    "binder": binder,
                    "binder_length": len(binder),
                    "pseudo_perplexity": ppl,
                    "has_unknown_aa": self.has_invalid_amino_acids(binder),
                }
            )
        return results

    def _sample_one(self, protein_seq: str, peptide_length: int, top_k: int, torch) -> str:
        """采样一条肽链。

        老脚本的做法是：一次性把 N 个 ``<mask>`` 拼在蛋白后面，取模型在
        **所有** mask 位置上的 logits，然后按 Top-K 逐位置独立采样。
        这里沿用同样的采样语义（各位置独立），但把非法氨基酸从候选里
        剔除。

        注意老脚本用 ``tokenizer.decode(predicted_token_ids)`` 把 token
        id 还原成序列——多 token 氨基酸（如 ``W``、``M`` 在某些词表里
        是双 token）会被拼回正确字符，但如果 decode 结果里混进空格，
        必须去掉，否则肽链长度与记录的长度不一致。
        """
        tokenizer = self._tokenizer
        masked_peptide = "<mask>" * peptide_length
        inputs = tokenizer(protein_seq + masked_peptide, return_tensors="pt").to(self._device)

        with torch.no_grad():
            logits = self._model(**inputs).logits

        mask_indices = (inputs["input_ids"] == tokenizer.mask_token_id).nonzero(as_tuple=True)[1]
        if mask_indices.numel() == 0:
            raise RuntimeError(
                "模型输出里找不到 <mask> 位置——tokenizer 可能不支持 <mask>，"
                "或输入序列被截断到了 0 长度。"
            )

        logits_at_masks = logits[0, mask_indices]
        top_k = min(top_k, logits_at_masks.shape[-1])

        top_k_logits, top_k_indices = logits_at_masks.topk(top_k, dim=-1)
        probabilities = torch.nn.functional.softmax(top_k_logits, dim=-1)

        if self.config.exclude_invalid_amino_acids and self._valid_token_ids is not None:
            top_k_indices, probabilities = self._filter_candidates(top_k_indices, probabilities, torch)

        # 用 torch 的采样保持与模型推理一致的分布
        predicted_indices = torch.multinomial(probabilities, num_samples=1)
        predicted_token_ids = top_k_indices.gather(-1, predicted_indices).squeeze(-1)

        binder = tokenizer.decode(predicted_token_ids, skip_special_tokens=True).replace(" ", "")
        return binder

    def _filter_candidates(self, top_k_indices, probabilities, torch):
        """把非法氨基酸的候选从 Top-K 里剔除并重新归一化概率。

        做法：构造一个布尔掩码，标记每个候选 token 是否是合法氨基酸，
        把非法候选的概率置零，然后重新归一化。这样采样永远不会选中
        ``X`` 之类的 token。

        如果某个位置上的候选**全部**非法（极端情况），则该位置退化为
        均匀分布——此时记录警告，而不是让 ``multinomial`` 抛错。
        """
        valid_ids = self._valid_token_ids or set()
        mask = torch.zeros_like(top_k_indices, dtype=torch.bool)
        for token_id in valid_ids:
            mask |= top_k_indices == token_id

        filtered_probs = probabilities * mask.to(probabilities.dtype)
        row_sums = filtered_probs.sum(dim=-1, keepdim=True)
        dead_rows = (row_sums.squeeze(-1) == 0)

        if bool(dead_rows.any()):
            LOGGER.warning(
                "有 %d 个位置的 Top-K 候选全部是非法氨基酸，已退化为均匀分布",
                int(dead_rows.sum().item()),
            )
            filtered_probs = torch.where(
                dead_rows.unsqueeze(-1), torch.ones_like(filtered_probs), filtered_probs
            )

        filtered_probs = filtered_probs / filtered_probs.sum(dim=-1, keepdim=True)
        return top_k_indices, filtered_probs

    # ------------------------------------------------------------------
    # 伪困惑度
    # ------------------------------------------------------------------
    def compute_pseudo_perplexity(self, protein_seq: str, binder_seq: str) -> float:
        """计算伪困惑度。

        做法：把「蛋白 + 肽链」一起编码，然后为肽链上的**每一个**位置
        分别做一次「只 mask 掉这一个位置」的推理，收集模型在该位置的
        交叉熵损失，取平均后取指数。

        老脚本用 ``masked_inputs = original_input.repeat(length_of_binder, 1)``
        一次性构造所有位置的掩码版本，效率更高。这里的实现保持同样的
        语义，但用 ``torch.nn.functional.cross_entropy`` 显式计算，
        避免依赖模型内部 ``labels`` 的平均方式（不同 transformers
        版本对 ``labels`` 的 reduction 行为并不完全一致）。

        Args:
            protein_seq: 靶标蛋白序列。
            binder_seq: 生成的肽链序列。

        Returns:
            伪困惑度值；越大表示模型越"意外"。
        """
        self._ensure_model()
        import torch

        tokenizer = self._tokenizer

        sequence = protein_seq + binder_seq
        input_ids = tokenizer.encode(sequence, return_tensors="pt").to(self._device)
        binder_length = len(binder_seq)

        if binder_length == 0:
            return float("inf")

        # 肽链对应的 token 区间（末尾 binder_length 个）
        positions = torch.arange(-binder_length, 0, device=self._device)

        losses: list[float] = []
        with torch.no_grad():
            for offset in range(binder_length):
                masked = input_ids.clone()
                masked[0, positions[offset]] = tokenizer.mask_token_id

                outputs = self._model(masked)
                logits = outputs.logits[0, positions[offset]]

                target = input_ids[0, positions[offset]]
                loss = torch.nn.functional.cross_entropy(
                    logits.unsqueeze(0), target.unsqueeze(0)
                )
                losses.append(float(loss.item()))

        mean_loss = sum(losses) / len(losses)
        return math.exp(mean_loss)

    # ------------------------------------------------------------------
    # 工具方法
    # ------------------------------------------------------------------
    @staticmethod
    def has_invalid_amino_acids(sequence: str) -> bool:
        """判断序列是否含有非标准氨基酸。"""
        return any(char not in VALID_AMINO_ACIDS for char in sequence)


def parse_target_records(
    fasta_path: Optional[Path] = None,
    csv_path: Optional[Path] = None,
) -> list[dict]:
    """读取靶标蛋白序列。

    支持 FASTA 与 CSV 两种输入。FASTA 的标识符提取规则沿用老脚本：
    优先取 ``GN=`` 字段（UniProt 的基因名），否则取 ``|`` 分隔的第二段
    （UniProt accession），否则取整行标题。

    Args:
        fasta_path: FASTA 文件路径。
        csv_path: CSV 文件路径，必须含 ``id`` 与 ``seq`` 两列。

    Returns:
        ``{"id": 标识符, "seq": 序列}`` 字典列表。
    """
    records: list[dict] = []
    if fasta_path:
        records.extend(_parse_fasta(fasta_path))
    if csv_path:
        records.extend(_parse_csv(csv_path))
    return records


def _parse_fasta(path: Path) -> list[dict]:
    """解析 FASTA，优先使用 ``GN=`` 字段作为标识符。"""
    records: list[dict] = []
    header: Optional[str] = None
    chunks: list[str] = []

    def flush() -> None:
        if header is None:
            return
        records.append({"id": _extract_id(header), "seq": "".join(chunks)})

    with open(path, encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                flush()
                header = line[1:]
                chunks = []
            elif header is not None:
                chunks.append(line.replace(" ", ""))
    flush()
    return records


def _extract_id(header: str) -> str:
    """从 FASTA 标题里提取标识符。

    优先 ``GN=``（基因名），其次 ``|`` 分隔的第二段（UniProt accession），
    最后退回标题首段。
    """
    if " GN=" in header:
        start = header.find(" GN=") + 4
        end = header.find(" ", start)
        if end == -1:
            end = header.find(" PE=", start)
        if end == -1:
            end = len(header)
        gn = header[start:end].strip()
        if gn:
            return gn

    parts = header.split("|")
    if len(parts) >= 2 and parts[1].strip():
        return parts[1].strip()

    return header.split()[0] if header.split() else header.strip()


def _parse_csv(path: Path) -> list[dict]:
    """解析 CSV 靶标文件，必须含 ``id`` 与 ``seq`` 列。"""
    import csv

    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return []
        columns = {name.strip().lower() for name in reader.fieldnames}
        if "id" not in columns or "seq" not in columns:
            raise ValueError(
                "CSV 靶标文件必须包含 'id' 与 'seq' 两列"
                f"（当前列: {', '.join(reader.fieldnames)}）"
            )
        id_key = next(name for name in reader.fieldnames if name.strip().lower() == "id")
        seq_key = next(name for name in reader.fieldnames if name.strip().lower() == "seq")

        records = []
        for row in reader:
            seq = (row.get(seq_key) or "").replace(" ", "").replace("\n", "").strip()
            if not seq:
                continue
            records.append({"id": (row.get(id_key) or "").strip(), "seq": seq})
        return records


def sample_from_probabilities(
    probabilities: Sequence[float],
    rng: Optional[random.Random] = None,
) -> int:
    """按概率分布采样一个下标（纯 Python 实现，用于测试与降级）。

    Args:
        probabilities: 概率序列，不必归一化。
        rng: 随机数发生器；None 时使用模块级 ``random``。

    Returns:
        采样到的下标。
    """
    generator = rng or random
    total = float(sum(probabilities))
    if total <= 0:
        return generator.randrange(len(probabilities))
    threshold = generator.random() * total
    cumulative = 0.0
    for index, probability in enumerate(probabilities):
        cumulative += float(probability)
        if cumulative >= threshold:
            return index
    return len(probabilities) - 1
