# skills-lock.json

本文件由 [`skills`](https://github.com/vercel-labs/skills) CLI 生成，记录
**这个仓库自己的技能**在注册表中的安装状态与内容哈希。

## 它解决什么问题

Agent Skills 分发的默认路径是 `npx skills add <owner>/<repo>`，安装的是
**仓库的当前 HEAD**——也就是说，任何一次推送都会改变用户装到的内容。

这在科研工具里是真实风险：用户按文档里的参数跑脚本，某个技能被改了一行，
结果与论文里的不一致，而用户无从察觉。

`skills-lock.json` 把这件事固定下来：

- `commit` —— 安装的是哪个 commit
- `computedHash` —— 该技能目录内容的 SHA-256（对文件做 CRLF→LF 归一后计算，
  因此在 Windows 与 Linux 上得到相同的值）

想复现本文档记录的版本：

```bash
npx skills add ssinder/bioinformatics-agent-skills --skill bulk-rnaseq
# 或在 CI 中用冻结模式：
npx skills install --frozen-lockfile
```

## 更新方式

技能内容变更后重新生成：

```bash
npx skills add ssinder/bioinformatics-agent-skills --skill <name> --yes
```

## 注意

这个文件只覆盖**本仓库自己的技能**。用户安装别人的技能时，
锁文件由他们自己的项目维护——那是 `skills install --frozen-lockfile`
发挥作用的地方。

CI 里没有强制这个文件是最新的：它会在有人改动技能内容后短暂地过期，
而这不影响任何测试的正确性（测试直接读源码）。
