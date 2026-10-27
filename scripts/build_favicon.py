#!/usr/bin/env python3
"""从 SVG 生成多尺寸 favicon.ico。

为什么需要这个脚本
------------------
favicon 不是「设计一次就完」的东西：同一张图在 64px 下好看，
到 16px（浏览器标签页、书签栏）可能糊成一团。所以 favicon 必须
**按尺寸分别设计**，再打包进同一个 .ico。

本脚本做两件事：
1. 用 cairosvg 把 docs/favicon.svg 渲染成各尺寸 PNG；
2. **手工组装 ICO 二进制** —— 因为 Pillow 的 ICO 插件不支持
   ``append_images``（传入也只会写入第一张图，实测只得到 16×16）。

ICO 容器格式::

    ICONDIR      : reserved(2)=0, type(2)=1, count(2)=N
    ICONDIRENTRY : width(1), height(1), colors(1), reserved(1),
                   planes(2), bpp(2), bytesInRes(4), imageOffset(4)
    image data   : 各张图的原始字节（Vista+ 支持直接内嵌 PNG）

其中 width/height 为 0 表示 256px（字段只有 1 字节，装不下 256）。

用法::

    python scripts/build_favicon.py
"""

from __future__ import annotations

import struct
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
FAVICON_SVG = REPO_ROOT / "docs" / "favicon.svg"
FAVICON_ICO = REPO_ROOT / "docs" / "favicon.ico"

#: 需要打包进 ICO 的尺寸（覆盖标签页、书签栏、桌面快捷方式）。
SIZES = (16, 32, 48, 64, 128, 256)


def _render_pngs(svg_path: Path, sizes: tuple[int, ...]) -> list[tuple[int, bytes]]:
    """把 SVG 渲染成各尺寸的 PNG 字节。"""
    try:
        import cairosvg
    except ImportError as exc:  # pragma: no cover - 依赖缺失时的提示
        raise SystemExit(
            "需要 cairosvg：pip install cairosvg\n"
            "（cairosvg 依赖系统 cairo 库；Windows 可用 conda install cairo）"
        ) from exc

    rendered: list[tuple[int, bytes]] = []
    for size in sizes:
        png = cairosvg.svg2png(url=str(svg_path), output_width=size, output_height=size)
        rendered.append((size, png))
    return rendered


def build_ico(images: list[tuple[int, bytes]]) -> bytes:
    """把多张 PNG 组装成 ICO 二进制。"""
    entries = b""
    data = b""
    offset = 6 + 16 * len(images)

    for size, blob in images:
        dimension = 0 if size >= 256 else size  # 1 字节字段，0 表示 256
        entries += struct.pack(
            "<BBBBHHII",
            dimension,  # width
            dimension,  # height
            0,  # 调色板颜色数（真彩色写 0）
            0,  # 保留
            1,  # color planes
            32,  # bits per pixel
            len(blob),  # 图像字节数
            offset,  # 在文件中的偏移
        )
        data += blob
        offset += len(blob)

    header = struct.pack("<HHH", 0, 1, len(images))  # reserved, type=1(icon), count
    return header + entries + data


def main() -> int:
    if not FAVICON_SVG.exists():
        print(f"缺少源文件：{FAVICON_SVG}")
        return 1

    images = _render_pngs(FAVICON_SVG, SIZES)
    ico = build_ico(images)
    FAVICON_ICO.write_bytes(ico)

    print(f"已生成 {FAVICON_ICO.relative_to(REPO_ROOT)}")
    print(f"  尺寸：{', '.join(f'{s}×{s}' for s in SIZES)}")
    print(f"  大小：{len(ico) / 1024:.1f} KB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
