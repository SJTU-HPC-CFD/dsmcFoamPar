# Codex 历史会话修复指南

## 适用症状

- `codex resume --all` 显示 `No sessions yet`
- VS Code 的 Codex/ChatGPT 面板看不到本地历史
- `state_5.sqlite` 里能查到线程，但 picker 和面板仍为空
- 同一台机器上同时存在 Linux 侧和 Windows 侧 `.codex`，状态相互漂移

## 这次验证过的根因

### 1. 只改 `state_5.sqlite` 不够

会话启动时会从 `sessions/.../rollout-*.jsonl` 第一行的 `session_meta` 反向修正数据库。

如果只改 `threads` 表，不改对应 `rollout` 文件里的：

- `cwd`
- `source`
- `model_provider`

下次启动后数据库会被写回旧值。

### 2. picker 识别的是一套一致的“本地会话”字段

这次最终证明，Linux 侧 `codex resume --all` 能稳定列出历史的前提是：

- `cwd = /home/superxcx/.codex`
- `source = cli`
- `model_provider = OpenAI`
- `archived = 0`
- `has_user_event = 1`

并且 `rollout_path` 必须指向真实存在的 Linux 路径。

### 3. Linux 侧不要把 `cwd` 写成 `\\wsl$...`

在 Linux 的 `thread/list` 路径里，`\\wsl$...` 会被当成相对路径拼到 `codexHome` 后面，导致过滤异常。

Linux 侧最终应写成：

```text
/home/superxcx/.codex
```

### 4. `thread/list` 是判断修复是否成功的关键

不要只看面板或 picker。

先确认后端 `thread/list` 已经能返回会话，再去看 `codex resume --all` 或 VS Code 面板。

### 5. VS Code 面板和 CLI 可能不是同一套 `.codex`

这次机器上同时存在：

- Linux: `/home/superxcx/.codex`
- Windows: `/mnt/c/Users/13637/.codex`

必须先确认真正被谁在读，再决定修哪一边。

## 可行修复顺序

### 第一步：整目录备份

不要只备份单个 sqlite。

至少备份：

- `state_5.sqlite`
- `session_index.jsonl`
- `history.jsonl`
- `sessions/`
- `archived_sessions/`

### 第二步：确认目标会话文件存在

检查目标会话对应的：

```text
sessions/YYYY/MM/DD/rollout-....jsonl
```

如果只在 `archived_sessions/` 里，需要先复制回 `sessions/`。

### 第三步：同时修两层元数据

#### A. 修 `rollout-*.jsonl` 第一行 `session_meta`

至少统一这三个字段：

- `cwd`
- `source`
- `model_provider`

Linux 侧推荐写成：

```json
"cwd": "/home/superxcx/.codex"
"source": "cli"
"model_provider": "OpenAI"
```

#### B. 修 `state_5.sqlite` 的 `threads` 行

同步统一：

- `cwd`
- `source`
- `model_provider`
- `archived`
- `archived_at`
- `has_user_event`
- `rollout_path`

### 第四步：重建平面索引

重建：

- `session_index.jsonl`
- `history.jsonl`

它们的内容必须和 `threads` 表一致。

### 第五步：做一次 WAL checkpoint

避免记录只停在 WAL 里：

```bash
sqlite3 ~/.codex/state_5.sqlite 'PRAGMA wal_checkpoint(TRUNCATE);'
```

### 第六步：先验证后端 `thread/list`

如果 `thread/list` 返回空，说明还没修到位。

如果 `thread/list` 能返回线程，而 picker 仍为空，才说明问题在 picker 层。

### 第七步：最后验证 picker

```bash
codex resume --all
```

## 单条会话最小修复模板

适合先验证一条关键会话是否能恢复：

1. 找到对应 `rollout-*.jsonl`
2. 改第一行 `session_meta.payload.cwd/source/model_provider`
3. 改 `state_5.sqlite` 中对应 `threads` 行
4. 重建 `session_index.jsonl` 和 `history.jsonl`
5. `PRAGMA wal_checkpoint(TRUNCATE)`
6. 先测 `thread/list`
7. 再测 `codex resume --all`

## 批量修复模板

当单条验证成功后，再批量处理：

1. 从已知正确备份库中读取所有线程元数据
2. 批量复制 `sessions/`
3. 批量改每个 `rollout` 文件的 `session_meta`
4. 用统一规则重建当前 `threads` 表
5. 统一重建 `session_index.jsonl` 和 `history.jsonl`
6. 再做 checkpoint
7. 再测 `thread/list`
8. 最后测 `codex resume --all`

## 这次实际验证有效的判断标准

不是“sqlite 里有记录”，而是下面三层同时通过：

1. `state_5.sqlite` 里能查到线程
2. `thread/list` 能返回线程
3. `codex resume --all` picker 能列出来

只有三层都通过，才算历史真正修复完成。
