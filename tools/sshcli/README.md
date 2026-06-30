# SnowFold SSH CLI - GPU 服务器连接管理工具

## 设计目标
解决 SSH 终端不连续问题：自动重连、会话保持、快速切换多台服务器。

## 核心功能
1. **一键连接**：记住服务器配置（IP、端口、用户、密钥）
2. **会话保持**：tmux 后台运行，断网不丢
3. **多服务器管理**：同时管理 5090/4090/3090 等多台
4. **文件传输**：scp/rsync 快捷命令
5. **状态监控**：GPU 使用、磁盘、显存
6. **便携**：纯 Python，单文件无依赖

## 调研结论
- tmux 是工业标准 (Linux/macOS 原生)
- zellij 是 Rust 现代化替代 (跨平台好)
- Windows 端用 Windows Terminal 集成
- VSCode 用户用 Remote-SSH 扩展

**我们的方案**: Python CLI 包装 tmux (无需额外学习) + 可选 zellij 后端
