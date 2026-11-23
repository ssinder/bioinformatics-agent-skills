# 开发常用命令
#
# 首次使用：make setup
# 提交前：make check

PYTHON ?= python
PIP    ?= $(PYTHON) -m pip

.DEFAULT_GOAL := help
.PHONY: help setup test lint validate check docs clean build release-note

help:  ## 显示所有可用命令
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-14s\033[0m %s\n", $$1, $$2}'

setup:  ## 安装开发依赖
	$(PIP) install -e ".[dev]"
	@echo "环境就绪。运行 'make check' 确认一切正常。"

test:  ## 运行测试套件
	$(PYTHON) -m pytest tests -q

lint:  ## 静态检查
	$(PYTHON) -m compileall -q skills packages scripts
	-$(PYTHON) -m ruff check skills packages scripts tests

validate:  ## 校验技能规范与仓库一致性
	$(PYTHON) scripts/validate_skills.py --strict
	$(PYTHON) scripts/ci_version_sync.py
	$(PYTHON) scripts/ci_consistency_check.py
	$(PYTHON) scripts/ci_check_fixtures.py

docs:  ## 重新生成文档索引
	$(PYTHON) scripts/build_docs_index.py

check: validate lint test  ## 提交前完整检查（规范 + 静态 + 测试）

build:  ## 构建 wheel 与 sdist
	$(PYTHON) -m build

clean:  ## 清理构建产物
	rm -rf build dist *.egg-info .pytest_cache
	find . -type d -name __pycache__ -not -path "./.git/*" -exec rm -rf {} +

release-note:  ## 打印本次发布相对上一个 tag 的变更摘要
	@$(PYTHON) -c "import subprocess; \
	last=$$(git describe --tags --abbrev=0 2>/dev/null || echo ''); \
	print(subprocess.run(['git','log','--oneline',f'{last}..HEAD'],capture_output=True,text=True).stdout)"
