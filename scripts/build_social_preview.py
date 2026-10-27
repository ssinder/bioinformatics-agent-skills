#!/usr/bin/env python3
"""从 SVG 源文件渲染社交预览图的 PNG 版本。

社交预览图（GitHub Social preview）只接受位图格式，SVG 不能直接上传。
本脚本负责把 ``docs/assets/social-preview.svg`` 渲染成 PNG，并同步生成
仓库根下的 ``.github/social-preview.png``（供外部工具直接取用）。

为什么要有这个脚本
------------------
上一版预览图是用一条临时命令渲染的，图上写死了当时的测试数
（829），测试数增长后就过时了，而且没有留下任何重新生成的手段 ——
想改个数字都得从头考古。把它固化成脚本后，改数据只需要改 SVG 源文件，
然后重跑本脚本。

用法
----
    python scripts/build_social_preview.py

依赖
----
    cairosvg（渲染）与 Pillow（PNG 校验与尺寸断言）

在 Windows 上若 pip 装到了错误的解释器，用 ``python -m pip install`` 确认
装到了当前 python 对应的环境里 —— 两条命令的 python 可能不是同一个。
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SVG_SOURCE = REPO_ROOT / "docs" / "assets" / "social-preview.svg"
PNG_TARGETS = (
    REPO_ROOT / "docs" / "assets" / "social-preview.png",
    REPO_ROOT / ".github" / "social-preview.png",
)

#: GitHub 社交预览图会按 2:1 裁切显示，宽高比不是 2:1 的图会被切掉两端。
EXPECTED_SIZE = (1280, 640)


def render(svg_path: Path, png_path: Path) -> None:
    """把 SVG 渲染成 PNG。"""
    import cairosvg

    cairosvg.svg2png(
        url=str(svg_path),
        write_to=str(png_path),
        output_width=EXPECTED_SIZE[0],
        output_height=EXPECTED_SIZE[1],
    )


def verify(png_path: Path) -> None:
    """校验渲染结果的尺寸 —— 比例不对会被 GitHub 裁切。"""
    from PIL import Image

    with Image.open(png_path) as im:
        actual = im.size
        mode = im.mode
    if actual != EXPECTED_SIZE:
        raise SystemExit(
            f"{png_path.name} 尺寸为 {actual}，期望 {EXPECTED_SIZE}；"
            "比例不符会导致 GitHub 裁切掉两端"
        )
    print(f"  [OK] {png_path.relative_to(REPO_ROOT)}  {actual[0]}x{actual[1]} {mode}")


def main() -> int:
    if not SVG_SOURCE.exists():
        raise SystemExit(f"找不到 SVG 源文件：{SVG_SOURCE}")

    print(f"从 {SVG_SOURCE.relative_to(REPO_ROOT)} 渲染社交预览图：")
    for target in PNG_TARGETS:
        target.parent.mkdir(parents=True, exist_ok=True)
        render(SVG_SOURCE, target)
        verify(target)

    print(f"\n完成，共生成 {len(PNG_TARGETS)} 个 PNG。")
    print("提示：GitHub Social preview 只接受位图，请上传 PNG 版本。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
