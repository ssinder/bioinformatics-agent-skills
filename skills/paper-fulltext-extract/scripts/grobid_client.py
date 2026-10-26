"""GROBID 全文抽取服务的 HTTP 客户端。

GROBID 是一个本地部署的文献结构化服务（Java），提供
``/api/processFulltextDocument`` 接口，把 PDF 转成 TEI XML。

**它不是云服务** —— 需要自己在本机或服务器上部署 GROBID（通常是
``docker run -it -p 8070:8070 kermitt/grobid``）。本模块只负责
HTTP 交互，不负责部署。

修掉的老仓库问题
----------------
原 ``paper_extract/scripts/batch_process_pdf.py`` 里有 5 个真实缺陷：

1. **URL 拼接脆弱** —— 默认值是 ``http://127.0.0.1:8070/``（带尾斜杠），
   代码直接拼 ``f"{url}api/processFulltextDocument"``。用户若传
   ``http://host:8070``（无尾斜杠），拼出来就是
   ``http://host:8070api/processFulltextDocument`` —— 404，且报错信息
   完全指不到真正原因。现在用 :func:`join_url` 统一处理。
2. **5xx 不重试** —— 只有超时才重试。GROBID 在处理大批量 PDF 时
   经常返回 503（过载）或 500，这类错误恰恰**应该**重试。
3. **重试耗尽后不计数** —— 三次超时后 `while` 循环自然结束，
   但既没 `fail_count += 1` 也没写失败记录，统计结果会撒谎。
4. **硬编码个人路径** —— ``--input`` / ``--output`` 的默认值是
   ``C:\\Users\\DELL\\Desktop\\trae\\paperdata\\...``。
5. **``action="store_true"`` 配 ``default=True``** —— 想关掉
   ``consolidateHeader`` 根本做不到。

本模块用零额外依赖实现（只用标准库 + urllib），requests 可选。
"""

from __future__ import annotations

import argparse
import json
import mimetypes
import os
import sys
import time
import urllib.error
import urllib.request
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional, Union

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from biopipe_core import cli_parser, get_logger, run_cli  # noqa: E402

LOGGER = get_logger("biopipe.paper-fulltext-extract")


__all__ = [
    "GrobidError",
    "GrobidResult",
    "GrobidClient",
    "join_url",
]


class GrobidError(RuntimeError):
    """GROBID 调用失败。"""


def join_url(base: str, path: str) -> str:
    """拼接 URL，正确处理尾斜杠。

    ``join_url("http://h:8070", "api/x")`` → ``http://h:8070/api/x``
    ``join_url("http://h:8070/", "api/x")`` → ``http://h:8070/api/x``

    这修掉了老代码里"默认带斜杠所以能用、用户传不带斜杠就 404"的问题。
    """
    return f"{base.rstrip('/')}/{path.lstrip('/')}"


@dataclass
class GrobidResult:
    """单篇文献的处理结果。"""

    pdf_path: str
    xml: Optional[str] = None
    elapsed_s: float = 0.0
    attempts: int = 0
    status: str = "pending"
    error: Optional[str] = None
    http_status: Optional[int] = None

    @property
    def ok(self) -> bool:
        return self.status == "ok" and bool(self.xml)


@dataclass
class GrobidClient:
    """GROBID HTTP 客户端。

    参数
    ----
    server
        服务地址，例如 ``http://127.0.0.1:8070``。
    timeout
        单次请求超时秒数。GROBID 处理一篇 PDF 通常 5-60 秒，
        复杂排版可能更久，默认给到 300 秒。
    max_retries
        失败重试次数（含首次尝试）。只对**可重试**的错误重试：
        网络异常、429、502、503、504。
    retry_delay
        首次重试前的等待秒数，之后按指数退避增长。
    consolidate_header
        合并参考文献中的作者/机构信息（GROBID 的
        ``consolidateHeader``）。标题质量问题多半出在这里没开。
    consolidate_citations
        合并同一文献的重复引用条目。
    flavor
        ``precision``（准，耗时更长）或 ``speed``（快，精度略低）。
    """

    server: str = "http://127.0.0.1:8070"
    timeout: float = 300.0
    max_retries: int = 3
    retry_delay: float = 5.0
    consolidate_header: bool = True
    consolidate_citations: bool = True
    flavor: str = "precision"
    extra_params: dict[str, str] = field(default_factory=dict)

    # 视为"可重试"的 HTTP 状态码
    RETRYABLE_STATUS = frozenset({429, 502, 503, 504})

    def __post_init__(self) -> None:
        self.server = self.server.rstrip("/")

    # ------------------------------------------------------------ 探测

    def is_available(self) -> tuple[bool, str]:
        """检查服务是否可达。

        返回 (可用, 说明)。调用 ``/api/isalive`` —— GROBID 提供的
        健康检查端点。
        """
        url = join_url(self.server, "api/isalive")
        try:
            with urllib.request.urlopen(url, timeout=10) as response:
                body = response.read().decode("utf-8", errors="replace").strip()
            if body.lower() == "true":
                return True, f"GROBID 服务可用（{self.server}）"
            return False, f"GROBID 服务返回异常状态: {body!r}"
        except urllib.error.HTTPError as exc:
            return False, f"GROBID 服务返回 HTTP {exc.code}"
        except urllib.error.URLError as exc:
            return False, f"无法连接 GROBID 服务 {self.server}: {exc.reason}"
        except Exception as exc:  # noqa: BLE001
            return False, f"检查 GROBID 服务时出错: {exc}"

    def version(self) -> str:
        """返回 GROBID 版本字符串；拿不到时返回 ``"unknown"``。"""
        url = join_url(self.server, "api/version")
        try:
            with urllib.request.urlopen(url, timeout=10) as response:
                return response.read().decode("utf-8", errors="replace").strip()
        except Exception:  # noqa: BLE001
            return "unknown"

    # ------------------------------------------------------------ 参数

    def build_params(self) -> dict[str, str]:
        """构造 ``processFulltextDocument`` 的表单参数。

        GROBID 期望的是**字符串**形式的布尔值（``"true"`` / ``"false"``），
        直接传 Python bool 会被编码成 ``"True"`` / ``"False"``，
        GROBID 解析时按非空即真处理，行为诡异但不报错 —— 很难查。
        """
        params = {
            "consolidateHeader": str(bool(self.consolidate_header)).lower(),
            "consolidateCitations": str(bool(self.consolidate_citations)).lower(),
            "flavor": self.flavor,
        }
        params.update(self.extra_params)
        return params

    # ------------------------------------------------------------ 请求

    def process_pdf(self, pdf_path: Union[str, os.PathLike]) -> GrobidResult:
        """把一个 PDF 交给 GROBID 处理，返回 TEI XML 文本。

        失败时**不抛异常**，而是把错误信息放进返回对象的 ``error`` 字段
        —— 批量处理时单篇失败不应中断整批。
        """
        path = Path(pdf_path)
        result = GrobidResult(pdf_path=str(path))

        if not path.exists():
            result.status = "skipped"
            result.error = f"文件不存在: {path}"
            return result

        params = self.build_params()
        url = join_url(self.server, "api/processFulltextDocument")
        boundary = uuid.uuid4().hex
        body = self._encode_multipart(path, params, boundary)

        start = time.monotonic()
        last_error = ""
        for attempt in range(1, max(1, self.max_retries) + 1):
            result.attempts = attempt
            try:
                request = urllib.request.Request(
                    url,
                    data=body,
                    headers={
                        "Content-Type": f"multipart/form-data; boundary={boundary}",
                        "Accept": "application/xml",
                    },
                    method="POST",
                )
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    payload = response.read().decode("utf-8", errors="replace")
                result.xml = payload
                result.status = "ok"
                result.http_status = 200
                result.elapsed_s = time.monotonic() - start
                return result

            except urllib.error.HTTPError as exc:
                result.http_status = exc.code
                detail = ""
                try:
                    detail = exc.read().decode("utf-8", errors="replace")[:300]
                except Exception:  # noqa: BLE001
                    pass
                last_error = f"HTTP {exc.code}: {detail}"

                if exc.code not in self.RETRYABLE_STATUS:
                    # 4xx 之类的客户端错误重试也没意义
                    break

            except urllib.error.URLError as exc:
                last_error = f"网络错误: {exc.reason}"
            except Exception as exc:  # noqa: BLE001
                last_error = f"{type(exc).__name__}: {exc}"
                break

            if attempt < self.max_retries:
                delay = self.retry_delay * (2 ** (attempt - 1))
                time.sleep(delay)

        result.status = "failed"
        result.error = last_error or "未知错误"
        result.elapsed_s = time.monotonic() - start
        return result

    def _encode_multipart(
        self,
        pdf_path: Path,
        params: dict[str, str],
        boundary: str,
    ) -> bytes:
        """把 PDF 与表单参数编码成 multipart/form-data。

        GROBID 要求文件字段名为 ``input``。
        """
        ctype = mimetypes.guess_type(pdf_path.name)[0] or "application/pdf"
        chunks: list[bytes] = []

        for key, value in params.items():
            chunks.append(
                f'--{boundary}\r\n'
                f'Content-Disposition: form-data; name="{key}"\r\n\r\n'
                f"{value}\r\n".encode()
            )

        chunks.append(
            f'--{boundary}\r\n'
            f'Content-Disposition: form-data; name="input"; '
            f'filename="{pdf_path.name}"\r\n'
            f"Content-Type: {ctype}\r\n\r\n".encode()
        )
        chunks.append(pdf_path.read_bytes())
        chunks.append(f"\r\n--{boundary}--\r\n".encode())
        return b"".join(chunks)

    # ------------------------------------------------------------ 批量

    def process_directory(
        self,
        input_dir: Union[str, os.PathLike],
        output_dir: Union[str, os.PathLike],
        skip_existing: bool = True,
        pause: float = 0.0,
        on_progress: Optional[Any] = None,
    ) -> list[GrobidResult]:
        """批量处理目录下的所有 PDF。

        参数
        ----
        input_dir
            PDF 所在目录（**非递归**）。
        output_dir
            TEI XML 输出目录。
        skip_existing
            已有对应 ``<name>.tei.xml`` 时跳过（断点续跑）。
        pause
            每篇之间的等待秒数，避免把服务打满。
        on_progress
            进度回调，形如 ``f(result, index, total)``。
        """
        src = Path(input_dir)
        dst = Path(output_dir)
        if not src.is_dir():
            raise GrobidError(f"输入目录不存在: {src}")
        dst.mkdir(parents=True, exist_ok=True)

        pdfs = sorted(
            p for p in src.iterdir() if p.is_file() and p.suffix.lower() == ".pdf"
        )
        results: list[GrobidResult] = []

        for index, pdf in enumerate(pdfs, start=1):
            xml_path = dst / f"{pdf.stem}.tei.xml"
            if skip_existing and xml_path.exists():
                result = GrobidResult(
                    pdf_path=str(pdf),
                    xml=xml_path.read_text(encoding="utf-8", errors="replace"),
                    status="skipped",
                )
                results.append(result)
                if on_progress is not None:
                    on_progress(result, index, len(pdfs))
                continue

            result = self.process_pdf(pdf)
            if result.ok:
                xml_path.write_text(result.xml or "", encoding="utf-8")
            results.append(result)
            if on_progress is not None:
                on_progress(result, index, len(pdfs))
            if pause > 0 and index < len(pdfs):
                time.sleep(pause)

        return results

    def summary(self, results: list[GrobidResult]) -> dict[str, Any]:
        """汇总批量处理结果。"""
        counts: dict[str, int] = {}
        for item in results:
            counts[item.status] = counts.get(item.status, 0) + 1
        return {
            "total": len(results),
            "by_status": counts,
            "succeeded": sum(1 for r in results if r.ok),
            "failed": sum(1 for r in results if r.status == "failed"),
            "skipped": sum(1 for r in results if r.status == "skipped"),
            "elapsed_s": round(sum(r.elapsed_s for r in results), 2),
        }


def _build_parser() -> argparse.ArgumentParser:
    """构造命令行解析器。"""
    parser = cli_parser(
        prog="grobid-client",
        description="调用 GROBID 服务把 PDF 批量转成 TEI XML",
        skill="paper-fulltext-extract",
    )
    parser.add_argument("--input", metavar="DIR", help="PDF 输入目录")
    parser.add_argument("--output", metavar="DIR", help="TEI XML 输出目录")
    parser.add_argument("--server", default="http://127.0.0.1:8070",
                        help="GROBID 服务地址（默认 http://127.0.0.1:8070）")
    parser.add_argument("--timeout", type=float, default=300.0, help="单次请求超时秒数")
    parser.add_argument("--max-retries", type=int, default=3, help="失败重试次数")
    parser.add_argument("--retry-delay", type=float, default=5.0, help="首次重试等待秒数")
    parser.add_argument("--flavor", choices=("precision", "speed"), default="precision",
                        help="抽取精度：precision（准）/ speed（快）")
    parser.add_argument("--no-consolidate-header", action="store_true",
                        help="关闭 consolidateHeader")
    parser.add_argument("--no-consolidate-citations", action="store_true",
                        help="关闭 consolidateCitations")
    parser.add_argument("--pause", type=float, default=0.0,
                        help="每篇之间的等待秒数，避免打满服务")
    parser.add_argument("--check-only", action="store_true",
                        help="仅检查服务是否可用并退出")
    return parser


def _main(args: argparse.Namespace) -> int:
    """命令行入口：检查 GROBID 服务并批量处理 PDF。"""
    client = GrobidClient(
        server=args.server,
        timeout=args.timeout,
        max_retries=args.max_retries,
        retry_delay=args.retry_delay,
        flavor=args.flavor,
        consolidate_header=not args.no_consolidate_header,
        consolidate_citations=not args.no_consolidate_citations,
    )

    available, message = client.is_available()
    LOGGER.info("%s", message)
    if not available:
        LOGGER.error(
            "GROBID 服务不可用。请先部署：\n"
            "  docker run -it -p 8070:8070 kermitt/grobid\n"
            "  首次启动需要下载模型，耗时较久，请耐心等待日志出现 'Server started'"
        )
        return 1

    LOGGER.info("GROBID 版本: %s", client.version())

    if args.check_only:
        return 0

    if not args.input or not args.output:
        LOGGER.error("需要同时提供 --input 与 --output")
        return 1

    def progress(result: GrobidResult, index: int, total: int) -> None:
        marker = {"ok": "✓", "skipped": "-", "failed": "✗"}.get(result.status, "?")
        LOGGER.info("[%d/%d] %s %s", index, total, marker, Path(result.pdf_path).name)
        if result.error:
            LOGGER.warning("    %s", result.error)

    results = client.process_directory(
        args.input, args.output, pause=args.pause, on_progress=progress
    )
    stats = client.summary(results)
    LOGGER.info(
        "完成：共 %d 篇，成功 %d，跳过 %d，失败 %d，耗时 %.1f 秒",
        stats["total"], stats["succeeded"], stats["skipped"],
        stats["failed"], stats["elapsed_s"],
    )
    return 0 if stats["failed"] == 0 else 1


if __name__ == "__main__":
    sys.exit(run_cli(_main, parser=_build_parser(), skill="paper-fulltext-extract", version="1.0.0"))
