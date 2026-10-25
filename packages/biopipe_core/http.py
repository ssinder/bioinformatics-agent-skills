"""HTTP 会话与重试工具。

科研脚本大量访问 PubMed、OpenAlex、UniProt、PubChem、AlphaFold 等公共 API，
这些服务普遍有速率限制和偶发 5xx。本模块提供统一的：

- 指数退避重试
- 限速节流
- 随机 User-Agent（降低被判定为爬虫的概率）
- 统一的错误类型 :class:`HttpError`
"""

from __future__ import annotations

import logging
import random
import time
from typing import Any, Optional

__all__ = ["http_session", "HttpError", "get_json", "post_json"]

_LOG = logging.getLogger("biopipe.http")

DEFAULT_USER_AGENT = (
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/122.0 Safari/537.36 biopipe/1.0"
)

_RETRY_STATUS = {408, 425, 429, 500, 502, 503, 504}


class HttpError(RuntimeError):
    """HTTP 请求失败（已耗尽重试次数）。"""

    def __init__(self, message: str, status: Optional[int] = None, url: Optional[str] = None):
        super().__init__(message)
        self.status = status
        self.url = url


def http_session(
    user_agent: Optional[str] = None,
    min_interval: float = 0.34,
    timeout: int = 60,
    retries: int = 4,
    backoff: float = 2.0,
    extra_headers: Optional[dict[str, str]] = None,
):
    """创建一个带重试与限速的 ``requests.Session``。

    Args:
        user_agent: 自定义 User-Agent。
        min_interval: 两次请求之间的最小间隔秒数（0.34 约等于 NCBI 每秒 3 次）。
        timeout: 单次请求超时秒数。
        retries: 失败重试次数。
        backoff: 指数退避基数。
        extra_headers: 附加请求头。

    Returns:
        配置好的 ``requests.Session``（requests 未安装时回退为 None）。
    """
    try:
        import requests
    except ImportError:  # pragma: no cover - 依赖缺失时的降级
        _LOG.warning("未安装 requests，HTTP 功能不可用")
        return None

    session = requests.Session()
    session.headers.update(
        {
            "User-Agent": user_agent or DEFAULT_USER_AGENT,
            "Accept-Language": "en-US,en;q=0.9",
            **(extra_headers or {}),
        }
    )

    if min_interval > 0:
        session.request = _throttled(session.request, min_interval)  # type: ignore[method-assign]

    session.request = _retrying(session.request, retries, backoff)  # type: ignore[method-assign]
    return session


def _throttled(original, min_interval: float):
    """包装原 request 方法，加入最小请求间隔。"""

    def wrapper(*args, **kwargs):
        time.sleep(min_interval)
        return original(*args, **kwargs)

    return wrapper


def _retrying(original, retries: int, backoff: float):
    """包装原 request 方法，加入指数退避重试。"""

    def wrapper(*args, **kwargs):
        last_exc: Optional[Exception] = None
        for attempt in range(retries + 1):
            try:
                response = original(*args, **kwargs)
            except Exception as exc:  # 网络异常
                last_exc = exc
                response = None

            if response is not None and response.status_code not in _RETRY_STATUS:
                return response

            if attempt >= retries:
                if response is not None:
                    raise HttpError(
                        f"请求失败，HTTP {response.status_code}: {response.url}",
                        status=response.status_code,
                        url=response.url,
                    )
                raise HttpError(f"请求失败: {last_exc}") from last_exc

            delay = (backoff**attempt) + random.uniform(0, 0.5)
            _LOG.warning(
                "请求失败(第 %d 次)，%.1fs 后重试: %s",
                attempt + 1,
                delay,
                response.status_code if response is not None else last_exc,
            )
            time.sleep(delay)

        raise HttpError("请求失败: 未知错误")  # pragma: no cover

    return wrapper


def get_json(url: str, session=None, params: Optional[dict[str, Any]] = None, **kwargs) -> Any:
    """发起 GET 请求并解析 JSON。

    Args:
        url: 请求地址。
        session: 可复用的 session，为 None 时临时创建。
        params: 查询参数。
        **kwargs: 透传给 ``session.get``。

    Returns:
        解析后的 JSON 对象。
    """
    own = session is None
    session = session or http_session()
    if session is None:
        raise HttpError("requests 不可用，无法发起请求")
    try:
        response = session.get(url, params=params, timeout=kwargs.pop("timeout", 60), **kwargs)
        response.raise_for_status()
        return response.json()
    finally:
        if own:
            session.close()


def post_json(url: str, session=None, json_body: Any = None, **kwargs) -> Any:
    """发起 POST 请求并解析 JSON。

    Args:
        url: 请求地址。
        session: 可复用的 session。
        json_body: JSON 请求体。
        **kwargs: 透传给 ``session.post``。
    """
    own = session is None
    session = session or http_session()
    if session is None:
        raise HttpError("requests 不可用，无法发起请求")
    try:
        response = session.post(
            url, json=json_body, timeout=kwargs.pop("timeout", 60), **kwargs
        )
        response.raise_for_status()
        return response.json()
    finally:
        if own:
            session.close()
