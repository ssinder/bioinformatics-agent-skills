## 变更类型 / Type of change

- [ ] 新增技能 / New skill
- [ ] 改进现有技能 / Improve existing skill
- [ ] 修复 bug / Bug fix
- [ ] 基础设施（CI、共享工具层、文档）/ Infrastructure
- [ ] 文档 / Documentation

## 关联技能 / Related skill

<!-- 填技能名，如 bulk-rnaseq -->

## 变更内容 / What changed

<!-- 简要说明 -->

## 自检清单 / Checklist

提交前请确认已勾选全部项目：

- [ ] `python scripts/validate_skills.py --strict` 零错误零警告
- [ ] `python -m pytest tests -q` 全部通过
- [ ] 新增/修改的脚本能响应 `--help`
- [ ] 没有硬编码绝对路径
- [ ] 没有裸 `print()`（统一使用 `logger`）
- [ ] 若新增技能，已更新 `packages/biopipe/registry.py` 的 `REGISTRY`
- [ ] `references/methodology.md` 与 `references/troubleshooting.md` 已补齐（新增技能时）
- [ ] 算法是真实实现，不是占位符

## 科学正确性 / Scientific correctness

<!-- 如果本次改动涉及科学算法，请说明：
     1. 算法的科学依据
     2. 如何验证结果正确
     3. 与原实现相比有无行为变化 -->

## 相关文档 / Documentation

<!-- 是否需要同步更新 README / docs/ 目录 -->
