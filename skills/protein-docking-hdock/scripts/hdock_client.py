#!/usr/bin/env python3
"""HDOCK Server 的 HTTP 客户端。

对应老仓库 ``dock_hdock/scripts/`` 下的 6 个 selenium 脚本，以及
``crawl_database/scripts/`` 下的 ``07``–``13`` 号 HDOCK/HADDOCK 脚本。

**HDOCK 没有公开的 REST API**，只有网页界面。老脚本因此全部用 Selenium
驱动 Edge 浏览器：启动浏览器、定位元素、填表、点提交、等结果页、解析表格。

但实测发现，HDOCK 的提交页是一个**标准的 HTML form POST**——没有 CSRF
token、没有隐藏字段、没有 JS 生成的签名。因此完全可以用 ``requests``
直接提交，结果页的分数也在普通 HTML 表格里，可以直接解析。

本模块把 selenium 换成了 HTTP，把"提交"和"取结果"拆成两个独立操作
（HDOCK 服务端计算通常需要 5–30 分钟，把两件事绑在一个脚本里会导致
用户长时间看不到任何输出）。

从老仓库继承并修掉的问题
------------------------
1. **硬编码个人路径**。``C:\\Users\\DELL\\Desktop\\protein\\msedgedriver.exe``、
   ``C:\\Users\\DELL\\Desktop\\hdock\\protein website.txt``、
   ``/Users/jung/Desktop/...`` 遍布全部脚本。
2. **硬编码个人邮箱**。``lisenlikun@foxmail.com`` 被直接写进提交代码里，
   别人用这个技能会往作者的邮箱发通知。
3. **多线程并发提交**。老脚本用 ``threading.Thread`` 起 N 个浏览器实例
   同时提交，没有任何间隔。HDOCK 是公益服务器，密集提交既容易被限流，
   也会拖慢所有人的作业。
4. **结果解析靠 XPath**。``//th[text()='Docking Score']/following-sibling::td``
   这种选择器在页面改版时会静默失败，返回 ``None`` 而不报错。
5. **文件名带空格与中文**。``"by hdock and af3 score get picture.py"`` 这类
   文件名在跨平台 shell 里处处是坑。

HDOCK 提交协议
--------------
首页 ``http://hdock.phys.hust.edu.cn/`` 上有一个 ``<form>``，它的
``action`` 属性是空的——真正的目标地址由页面里的 JavaScript 在提交前
改写::

    function OnSubmitForm(){
        ...
        if( ! sym.value ){
            document.form1.action = "runhdock.php";     // 普通对接
        }else{
            document.form1.action = "runhsymdock.php";  // 对称多聚体对接
        }
        return true;
    }

**这是最容易踩的坑**：直接 POST 到首页会被当成普通页面请求，返回字节级
相同的首页 HTML，既不报错也没有结果。必须 POST 到 ``runhdock.php``。

表单字段（全部为 ``multipart/form-data``）::

    受体: pdbfile1 | pdbid1 | fastaseq1 | fastafile1   （四选一）
    配体: pdbfile2 | pdbid2 | fastaseq2 | fastafile2   （四选一）
    其它: email, jobname, ligtyp, symmetry,
          sitenum1, sitenum2, pdbsite, restrnum,
          sitefile1, sitefile2, restrfile, saxsfile

结果页有三种状态
----------------
1. **排队中** —— 页面含 ``is QUEUED``，服务端每 10 秒自动刷新一次。
2. **已完成** —— 页面含 ``<th>Docking Score</th>`` 与
   ``<th>Confidence Score</th>`` 及其兄弟 ``<td>``。
3. **失败** —— 页面含 ``ERROR:``，例如 ``No structure found for the
   Receptor molecule``。
"""

from __future__ import annotations

import re
import sys
import time
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional
from urllib.parse import urljoin

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))

from biopipe_core import get_logger, http_session  # noqa: E402

__all__ = [
    "HDOCK_BASE_URL",
    "HDOCK_SUBMIT_URL",
    "HDOCK_SYMMETRIC_SUBMIT_URL",
    "JobStatus",
    "DockingResult",
    "submit_docking_job",
    "fetch_result_page",
    "parse_result_page",
    "wait_for_results",
    "download_top_models",
    "compute_confidence_score",
    "LIGAND_TYPES",
    "build_form_fields",
]

LOGGER = get_logger("biopipe.protein-docking-hdock")

HDOCK_BASE_URL = "http://hdock.phys.hust.edu.cn/"
#: 普通对接的提交入口。首页 form 的 action 是空的，真实地址由 JS 改写。
HDOCK_SUBMIT_URL = urljoin(HDOCK_BASE_URL, "runhdock.php")
#: 对称多聚体对接（填了 symmetry 时）的提交入口。
HDOCK_SYMMETRIC_SUBMIT_URL = urljoin(HDOCK_BASE_URL, "runhsymdock.php")

#: 配体类型。取自首页 ``<select name="ligtyp">`` 的 option value。
LIGAND_TYPES = {
    "protein": "Protein",
    "rna": "ssRNA",
    "dna": "dsDNA",
    "rna.d": "dsRNA",
    "ssdna": "ssDNA",
}

#: 结果页里用 ``<th>``/``<td>`` 成对呈现的字段
_SCORE_PAIR_RE = re.compile(
    r"<th[^>]*>\s*(?P<label>[^<]+?)\s*</th>\s*<td[^>]*>\s*(?P<value>[^<]*?)\s*</td>",
    re.IGNORECASE | re.DOTALL,
)
#: 排名矩阵表的一整行：``<th>标签</th><td>值</td><td>值</td>...``
_SCORE_ROW_RE = re.compile(
    r"<tr[^>]*>\s*<th[^>]*>\s*(?P<label>[^<]+?)\s*</th>(?P<cells>(?:\s*<td[^>]*>.*?</td>)+)\s*</tr>",
    re.IGNORECASE | re.DOTALL,
)
_CELL_RE = re.compile(r"<td[^>]*>(.*?)</td>", re.IGNORECASE | re.DOTALL)
_ERROR_RE = re.compile(r"ERROR:\s*(.*?)(?:</h2>|</p>|<br|$)", re.IGNORECASE | re.DOTALL)
_JOB_URL_RE = re.compile(r"/data/([A-Za-z0-9]+)/")
_TAG_RE = re.compile(r"<[^>]+>")


class JobStatus:
    """作业状态常量。"""

    QUEUED = "queued"
    DONE = "done"
    FAILED = "failed"
    UNKNOWN = "unknown"


@dataclass
class DockingResult:
    """一个 HDOCK 作业的结果。"""

    job_name: str
    job_id: str = ""
    url: str = ""
    status: str = JobStatus.UNKNOWN
    docking_score: str = ""
    confidence_score: str = ""
    #: 结果页里所有 ``<th>``/``<td>`` 字段的原始键值对
    fields: dict[str, str] = field(default_factory=dict)
    #: 排名矩阵：``{"Docking Score": ["-444.95", "-244.21", ...]}``
    #: 列表按排名（1、2、3…）排列
    rankings: dict[str, list[str]] = field(default_factory=dict)
    errors: list[str] = field(default_factory=list)
    #: 可下载的结果文件链接（文件名 → 完整 URL）
    downloads: dict[str, str] = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        """作业是否成功完成。"""
        return self.status == JobStatus.DONE

    @property
    def confidence_score_computed(self) -> float:
        """按官方公式从 Docking Score 算出的 Confidence Score。"""
        if not self.docking_score:
            return 0.0
        return compute_confidence_score(self.docking_score)

    def to_row(self) -> dict:
        """转成扁平字典，供 CSV 导出。"""
        return {
            "job_name": self.job_name,
            "job_id": self.job_id,
            "url": self.url,
            "status": self.status,
            "docking_score": self.docking_score,
            "confidence_score": self.confidence_score,
            "error": "; ".join(self.errors),
        }


def build_form_fields(
    job_name: str = "",
    email: str = "",
    ligand_type: str = "protein",
    symmetry: str = "",
    receptor_pdb_id: str = "",
    ligand_pdb_id: str = "",
    receptor_seq: str = "",
    ligand_seq: str = "",
    receptor_site: str = "",
    ligand_site: str = "",
    restraint: str = "",
) -> dict[str, str]:
    """构造 HDOCK 表单字段。

    受体与配体各支持四种输入方式（PDB 文件 / PDB ID:ChainID / FASTA
    序列 / FASTA 文件），一次只能用一种。

    Args:
        job_name: 作业名（服务端会把它里的空白替换成下划线）。
        email: 通知邮箱。为空则服务端不会发通知。
        ligand_type: 配体类型，见 :data:`LIGAND_TYPES`。
        symmetry: 对称性（如 ``C2``、``D3``）。填了会走对称对接入口。
        receptor_pdb_id: 受体的 ``PDBID:ChainID``，如 ``1CGI:E``。
        ligand_pdb_id: 配体的 ``PDBID:ChainID``。
        receptor_seq: 受体 FASTA 序列。
        ligand_seq: 配体 FASTA 序列。
        receptor_site: 受体结合位点残基号，逗号分隔。
        ligand_site: 配体结合位点残基号，逗号分隔。
        restraint: 残基距离限制。

    Returns:
        可直接 POST 的字段字典。
    """
    fields = {
        "fastaseq1": receptor_seq,
        "fastaseq2": ligand_seq,
        "pdbid1": receptor_pdb_id,
        "pdbid2": ligand_pdb_id,
        "email": email,
        "jobname": job_name,
        "ligtyp": ligand_type,
        "symmetry": symmetry,
        "sitenum1": receptor_site,
        "sitenum2": ligand_site,
        "pdbsite": "",
        "restrnum": restraint,
        # 文件字段留空；上传文件时由调用方通过 files= 传入
        "pdbfile1": "",
        "pdbfile2": "",
        "fastafile1": "",
        "fastafile2": "",
        "sitefile1": "",
        "sitefile2": "",
        "restrfile": "",
        "saxsfile": "",
        "upload": "Submit",
    }
    return fields


def compute_confidence_score(docking_score: float) -> float:
    """按 HDOCK 官方公式从 Docking Score 计算 Confidence Score。

    HDOCK 的帮助文档给出了闭式公式::

        Confidence_score = 1.0 / [1.0 + e^(0.02 * (Docking_Score + 150))]

    这个分数是 docking score 的单调函数，所以它不提供 Docking Score
    之外的新信息——但**它提供了可直接比较的标量**：经验上 >0.7 表示
    两个分子很可能结合。

    结果页有时不给这一列，这时用这个函数补上，比在页面上等一个可能
    永远不出现的元素可靠得多。

    Args:
        docking_score: Docking Score（越负越好）。

    Returns:
        Confidence Score，取值范围 (0, 1)。
    """
    import math

    try:
        score = float(docking_score)
    except (TypeError, ValueError):
        return 0.0
    try:
        return 1.0 / (1.0 + math.exp(0.02 * (score + 150.0)))
    except OverflowError:  # pragma: no cover - score 极端负时
        return 1.0


def _file_tuple(path: Path) -> tuple:
    """构造 requests 的文件上传元组。"""
    handle = path.open("rb")
    return (handle, path.name)


def submit_docking_job(
    job_name: str = "",
    email: str = "",
    ligand_type: str = "protein",
    symmetry: str = "",
    receptor_pdb_id: str = "",
    ligand_pdb_id: str = "",
    receptor_seq: str = "",
    ligand_seq: str = "",
    receptor_pdb_file: Optional[str] = None,
    ligand_pdb_file: Optional[str] = None,
    receptor_fasta_file: Optional[str] = None,
    ligand_fasta_file: Optional[str] = None,
    receptor_site: str = "",
    ligand_site: str = "",
    restraint: str = "",
    session=None,
    timeout: int = 180,
) -> DockingResult:
    """向 HDOCK 提交一个对接作业。

    Args:
        job_name: 作业名。
        email: 通知邮箱。
        ligand_type: 配体类型，见 :data:`LIGAND_TYPES`。
        symmetry: 对称性；填了会改用对称对接入口。
        receptor_pdb_id: 受体 ``PDBID:ChainID``。
        ligand_pdb_id: 配体 ``PDBID:ChainID``。
        receptor_seq: 受体 FASTA 序列。
        ligand_seq: 配体 FASTA 序列。
        receptor_pdb_file: 受体 PDB 文件路径。
        ligand_pdb_file: 配体 PDB 文件路径。
        receptor_fasta_file: 受体 FASTA 文件路径。
        ligand_fasta_file: 配体 FASTA 文件路径。
        receptor_site: 受体结合位点残基号。
        ligand_site: 配体结合位点残基号。
        restraint: 残基距离限制。
        session: 可复用的 ``requests.Session``。
        timeout: 单次请求超时秒数（提交可能需要几十秒）。

    Returns:
        包含作业 URL 的 :class:`DockingResult`（``status`` 通常是 ``queued``）。
    """
    has_receptor = any(
        [receptor_pdb_id, receptor_seq, receptor_pdb_file, receptor_fasta_file]
    )
    has_ligand = any([ligand_pdb_id, ligand_seq, ligand_pdb_file, ligand_fasta_file])
    if not has_receptor:
        raise ValueError(
            "受体的输入缺失：需要提供 PDB 文件、PDB ID:ChainID、FASTA 序列"
            "或 FASTA 文件（四选一）"
        )
    if not has_ligand:
        raise ValueError(
            "配体的输入缺失：需要提供 PDB 文件、PDB ID:ChainID、FASTA 序列"
            "或 FASTA 文件（四选一）"
        )

    fields = build_form_fields(
        job_name=job_name,
        email=email,
        ligand_type=ligand_type,
        symmetry=symmetry,
        receptor_pdb_id=receptor_pdb_id,
        ligand_pdb_id=ligand_pdb_id,
        receptor_seq=receptor_seq,
        ligand_seq=ligand_seq,
        receptor_site=receptor_site,
        ligand_site=ligand_site,
        restraint=restraint,
    )

    files: dict[str, tuple] = {}
    for key, value in (
        ("pdbfile1", receptor_pdb_file),
        ("pdbfile2", ligand_pdb_file),
        ("fastafile1", receptor_fasta_file),
        ("fastafile2", ligand_fasta_file),
    ):
        if value:
            files[key] = _file_tuple(Path(value))

    own_session = session is None
    if own_session:
        session = http_session(min_interval=1.0, retries=3)
    if session is None:
        LOGGER.warning("requests 不可用，无法提交 HDOCK 作业")
        return DockingResult(job_name=job_name, status=JobStatus.FAILED)

    url = HDOCK_SYMMETRIC_SUBMIT_URL if symmetry else HDOCK_SUBMIT_URL

    try:
        # 先 GET 一次首页，建立与服务器的基本交互（HDOCK 不发 cookie，
        # 但这与浏览器行为一致，也确保 User-Agent 等头一致）
        session.get(HDOCK_BASE_URL, timeout=timeout)

        if files:
            response = session.post(
                url, data=fields, files=files, timeout=timeout, allow_redirects=True
            )
        else:
            response = session.post(
                url, data=fields, timeout=timeout, allow_redirects=True
            )
        response.raise_for_status()
    except Exception as exc:  # noqa: BLE001
        LOGGER.warning("HDOCK 提交失败: %s", exc)
        return DockingResult(
            job_name=job_name,
            status=JobStatus.FAILED,
            errors=[f"提交失败: {exc}"],
        )
    finally:
        if own_session:
            for handle, _name in files.values():
                handle.close()
            session.close()

    result_url = response.url
    match = _JOB_URL_RE.search(result_url)
    job_id = match.group(1) if match else ""

    LOGGER.info("作业已提交: %s", result_url)
    return DockingResult(
        job_name=job_name,
        job_id=job_id,
        url=result_url,
        status=JobStatus.QUEUED,
    )


def fetch_result_page(url: str, session=None, timeout: int = 60) -> str:
    """抓取结果页 HTML。

    Args:
        url: 作业结果页 URL。
        session: 可复用的 ``requests.Session``。
        timeout: 单次请求超时秒数。

    Returns:
        页面 HTML 文本；失败时返回空串。
    """
    own_session = session is None
    if own_session:
        session = http_session(min_interval=1.0, retries=3)
    if session is None:
        return ""

    try:
        response = session.get(url, timeout=timeout)
        response.raise_for_status()
        return response.text
    except Exception as exc:  # noqa: BLE001
        LOGGER.warning("获取结果页失败: %s", exc)
        return ""
    finally:
        if own_session:
            session.close()


def parse_result_page(html: str, job_name: str = "", url: str = "") -> DockingResult:
    """解析结果页。

    识别三种状态：排队中、已完成、失败。

    Args:
        html: 结果页 HTML。
        job_name: 作业名。
        url: 作业 URL。

    Returns:
        解析出的 :class:`DockingResult`。
    """
    result = DockingResult(job_name=job_name, url=url)
    if not html:
        result.status = JobStatus.UNKNOWN
        return result

    match = _JOB_URL_RE.search(url or "")
    if match:
        result.job_id = match.group(1)

    # 失败状态：ERROR: No structure found for the Receptor molecule.
    errors = [re.sub(r"<[^>]+>", "", e).strip() for e in _ERROR_RE.findall(html)]
    if errors:
        result.errors = [e for e in errors if e]
        result.status = JobStatus.FAILED
        return result

    # 排队状态
    if "is QUEUED" in html:
        result.status = JobStatus.QUEUED
        return result

    # 完成状态：结果页是一个**排名矩阵表**——行是指标，列是排名 1..N。
    # 例如：
    #   <tr><th>Docking Score</th><td>-444.95</td><td>-244.21</td>...</tr>
    #   <tr><th>Ligand rmsd (Å)</th><td>0.42</td><td>18.23</td>...</tr>
    # 取第一格作为该作业的"最佳"分数，同时保留整行用于排名分析。
    for label, cells_html in _SCORE_ROW_RE.findall(html):
        key = re.sub(r"\s+", " ", _TAG_RE.sub("", label)).strip()
        if not key:
            continue
        values = [
            re.sub(r"\s+", " ", _TAG_RE.sub("", cell)).strip()
            for cell in _CELL_RE.findall(cells_html)
        ]
        if not values:
            continue
        result.rankings[key] = values
        result.fields.setdefault(key, values[0])
        if key == "Docking Score":
            result.docking_score = values[0]
        elif key == "Confidence Score":
            result.confidence_score = values[0]

    # 有些结果页用非表格的形式给出 Confidence Score
    if not result.confidence_score:
        found = re.search(
            r"Confidence\s*Score[^0-9\-]{0,40}(-?\d+\.?\d*)", html, re.IGNORECASE
        )
        if found:
            result.confidence_score = found.group(1).strip()

    # 结果页有时不直接给 Confidence Score，但官方公式是公开的
    if result.docking_score and not result.confidence_score:
        computed = compute_confidence_score(result.docking_score)
        result.confidence_score = f"{computed:.4f}"
        result.fields.setdefault("Confidence Score (computed)", f"{computed:.4f}")

    if result.docking_score or result.confidence_score:
        result.status = JobStatus.DONE
    else:
        result.status = JobStatus.UNKNOWN

    # 下载链接：按 href 里的文件名识别（链接文本常常是 "Click to Show"）
    #
    # **注意 href 可能没有引号**。HDOCK 真实的结果页写的是
    # ``<a target=_blank href=model_1.pdb>``，而不是 ``href="model_1.pdb"``。
    # 只匹配带引号的形式会**一个模型都找不到**——而且不报错，只是
    # downloads 是个空字典。
    #
    # 这里用命名分组 + 非捕获组，让 findall 始终返回 (href, text) 两元组。
    link_re = re.compile(
        r'<a[^>]*?href=(?:"(?P<dq>[^"]+)"|\'(?P<sq>[^\']+)\'|(?P<bare>[^\s>"\']+))'
        r'[^>]*>(?P<text>.*?)</a>',
        re.IGNORECASE | re.DOTALL,
    )
    for match in link_re.finditer(html):
        href = match.group("dq") or match.group("sq") or match.group("bare") or ""
        if not href:
            continue
        clean = re.sub(r"\s+", " ", _TAG_RE.sub("", match.group("text"))).strip()
        # href 里的文件名才是真正的资源名
        name = href.rstrip("/").split("/")[-1]
        if not name:
            continue
        if re.search(r"\.(pdb|zip|out|txt|cif|json)$", name, re.IGNORECASE):
            label = clean if re.search(r"\.(pdb|zip|out|txt|cif|json)$", clean, re.IGNORECASE) else name
            result.downloads[label] = urljoin(url or HDOCK_BASE_URL, href)

    return result


def wait_for_results(
    url: str,
    job_name: str = "",
    poll_interval: int = 20,
    max_wait: int = 3600,
    session=None,
    timeout: int = 60,
    on_poll=None,
) -> DockingResult:
    """轮询等待作业完成。

    HDOCK 服务端计算通常需要 5–30 分钟，官方文档写的是"通常 30 分钟内
    完成"。这里默认最多等 1 小时。

    Args:
        url: 作业结果页 URL。
        job_name: 作业名。
        poll_interval: 轮询间隔秒数（服务端每 10 秒刷新一次页面，
            没必要比这更频繁）。
        max_wait: 最长等待秒数。
        session: 可复用的 ``requests.Session``。
        timeout: 单次请求超时秒数。
        on_poll: 每次轮询后的回调，接收当前的 :class:`DockingResult`。

    Returns:
        最终的 :class:`DockingResult`。超时返回状态为 ``queued`` 的结果。
    """
    deadline = time.time() + max_wait
    attempt = 0
    while time.time() < deadline:
        attempt += 1
        html = fetch_result_page(url, session=session, timeout=timeout)
        result = parse_result_page(html, job_name=job_name, url=url)
        if on_poll is not None:
            on_poll(result)
        if result.status in (JobStatus.DONE, JobStatus.FAILED):
            return result
        LOGGER.info(
            "作业仍在排队（第 %d 次检查，每 %d 秒一次）", attempt, poll_interval
        )
        time.sleep(poll_interval)

    LOGGER.warning("等待超时（%d 秒），作业可能仍在服务端计算", max_wait)
    return DockingResult(
        job_name=job_name,
        url=url,
        status=JobStatus.QUEUED,
        errors=[f"等待超过 {max_wait} 秒仍未完成"],
    )


def download_top_models(
    url: str,
    output_dir: Path,
    top_n: int = 5,
    session=None,
    timeout: int = 120,
) -> list[Path]:
    """下载排名靠前的复合物模型。

    Args:
        url: 作业结果页 URL。
        output_dir: 输出目录。
        top_n: 下载前 N 个模型。
        session: 可复用的 ``requests.Session``。
        timeout: 单次请求超时秒数。

    Returns:
        实际下载到的文件路径列表。
    """
    html = fetch_result_page(url, session=session, timeout=timeout)
    result = parse_result_page(html, url=url)
    if result.status != JobStatus.DONE:
        LOGGER.warning("作业尚未完成，无法下载模型")
        return []

    own_session = session is None
    if own_session:
        session = http_session(min_interval=1.0, retries=3)
    if session is None:
        return []

    output_dir.mkdir(parents=True, exist_ok=True)
    downloaded: list[Path] = []
    try:
        # 按文件名里的排名数字排序：top_1.pdb、top_2.pdb …
        def _rank(item: tuple[str, str]) -> int:
            found = re.findall(r"(\d+)", item[0])
            return int(found[0]) if found else 9999

        for name, link in sorted(result.downloads.items(), key=_rank):
            if len(downloaded) >= top_n:
                break
            if not name.lower().endswith(".pdb"):
                continue
            try:
                response = session.get(link, timeout=timeout)
                response.raise_for_status()
                target = output_dir / name
                target.write_bytes(response.content)
                downloaded.append(target)
                LOGGER.info("已下载: %s", target)
            except Exception as exc:  # noqa: BLE001
                LOGGER.warning("下载 %s 失败: %s", name, exc)
    finally:
        if own_session:
            session.close()

    return downloaded
