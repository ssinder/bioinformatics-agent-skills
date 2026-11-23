"""NCBI E-utilities 调用身份与免责声明。

为什么需要这个模块
------------------
NCBI 对 E-utilities 的使用有明确规定
（https://www.ncbi.nlm.nih.gov/books/NBK25497/），其中三条直接影响
本项目的设计：

1. **email 必须是软件开发者��邮箱，而不是终端用户的邮箱。**
   原文："The value of **email** should be a complete and valid e-mail
   address of the software developer and not that of a third-party end
   user."
   → 正确做法是由项目方注册并内置一个身份，而不是要求每个用户
   自己去导出环境变量。

2. **仅在请求里带上 tool/email 并不算合规**，这两个值必须事先向
   NCBI 注册。原文："merely providing values for `tool` and `email` in
   requests is not sufficient to comply with this policy; these values
   must be registered with NCBI."

3. **在软件中使用 E-utilities 时，NCBI 的免责声明与版权声明必须对
   用户可见。** 原文："If you use the E-utilities within software,
   NCBI's Disclaimer and Copyright notice must be evident to users of
   your product."

另外，限速规则是：无 API key 时超过 3 req/s 会报错，有 key 时可到
10 req/s。也就是说 **API key 是可选的性能开关，不是使用前提**。

身份配置的优先级
----------------
1. 命令行 ``--email``（用户显式覆盖）
2. 环境变量 ``BIOAGENT_EMAIL``
3. 本模块内置的**项目方注册身份**（默认）
4. 环境变量 ``NCBI_EMAIL``（兼容常见的通用命名）
"""

from __future__ import annotations

import os
from typing import Optional

__all__ = [
    "NCBI_TOOL",
    "DEFAULT_EMAIL",
    "NCBI_DISCLAIMER",
    "resolve_email",
    "disclaimer_block",
]

# NCBI 要求的软件标识：唯一标识产生请求的软件，不含内部空格。
NCBI_TOOL = "bioinformatics-agent-skills"

# 项目方注册的开发者联系邮箱。
#
# 说明：NCBI 要求这个地址属于**软件开发者**而非终端用户。维护者需要在
# 这里换成自己收到通知的邮箱，并在 NCBI 处完成注册——注册的流程见
# https://www.ncbi.nlm.nih.gov/books/NBK25497/ 。
#
# 之所以仍然允许用户覆盖，是因为部分用户所在机构的 IT 政策要求
# 出口流量使用本机构身份。
DEFAULT_EMAIL = "bioinformatics-agent-skills@users.noreply.github.com"

# NCBI 要求的免责声明与版权声明。必须在使用 E-utilities 的产品中对
# 用户可见——本模块把它接到了各技能脚本的输出头部。
NCBI_DISCLAIMER = (
    "NCBI Disclaimer and Copyright Notice:\n"
    "  This software uses the NCBI E-utilities. The data obtained through\n"
    "  the E-utilities are subject to the NCBI Disclaimer and Copyright\n"
    "  Notice: https://www.ncbi.nlm.nih.gov/about/disclaimer.html\n"
    "  The NCBI databases are made available for use by individual researchers\n"
    "  and for use in the design, development, and implementation of software\n"
    "  and databases. The data are subject to copyright protection by the\n"
    "  National Library of Medicine and the U.S. government. Portions of\n"
    "  the data may be subject to additional restrictions by their original\n"
    "  producers. For more information, see\n"
    "  https://www.ncbi.nlm.nih.gov/books/NBK25497/."
)


def resolve_email(explicit: Optional[str] = None) -> str:
    """解析本次请求使用的 NCBI 联系邮箱。

    优先级：显式参数 > ``BIOAGENT_EMAIL`` > 内置的项目方注册身份 >
    ``NCBI_EMAIL``。

    永远返回一个非空字符串——内置身份保证了这一点，所以用户**不需要**
    为了跑通脚本去做任何配置。
    """
    for candidate in (
        explicit,
        os.environ.get("BIOAGENT_EMAIL"),
        DEFAULT_EMAIL,
        os.environ.get("NCBI_EMAIL"),
    ):
        if candidate and candidate.strip():
            return candidate.strip()
    return DEFAULT_EMAIL


def disclaimer_block() -> str:
    """返回应展示给用户的 NCBI 免责声明文本。"""
    return NCBI_DISCLAIMER
