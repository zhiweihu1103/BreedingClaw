<div align="center">
<img src='https://github.com/user-attachments/assets/e62ee276-c4e9-4ca8-92ea-5779dff4ea22' alt='BreedingClaw' height='150px'>
</div>

# <p align="center"><font size=50><strong>[首个育种领域小龙虾](https://github.com/zhiweihu1103/BreedingClaw)</strong></font></p>

<div align="center"> <img src='https://github.com/user-attachments/assets/aa11d15e-f55a-4b48-af90-b8be6e101a87' alt='山西省后稷实验室' height='90px'> &emsp;&emsp;&emsp;&emsp; <!-- Adds extra spacing between images --> <img src='https://github.com/user-attachments/assets/4d0aa98e-8b86-4d34-b8cd-1f65ec6256ed' alt='山西农业大学' height='90px'> &emsp;&emsp;&emsp;&emsp; <!-- Adds extra spacing between images --> <img src='https://github.com/user-attachments/assets/94aeaa28-45dc-446c-9a8d-9f39fc51bd44' alt='山西大学' height='90px'> &emsp;&emsp;&emsp;&emsp; <!-- Adds extra spacing between images --> 
</div>

---
## 目录

- [概览](#概览)
- [快速开始](#快速开始)
- [消息通道](#消息通道)

## 概览

**BreedingClaw**是育种领域的第一只龙虾，专注于**育种基因的生物信息分析**。它通过自然语言交互，让研究者可以完成以下任务：

- 基因序列比对与 BLAST 检索  
- 育种候选基因分析与可视化  
- 蛋白质结构渲染（PyMOL）  
- 测序数据质控和差异分析（FastQC / MultiQC / 火山图等）  
- 文献检索与摘要

## 快速开始

**1.安装**

**Install from source**

```bash
git clone https://github.com/zhiweihu1103/BreedingClaw.git
cd BreedingClaw
pip install -e .
```
**Install from PyPI** 

```bash
pip install nanobot-ai
```

**2.初始化**

```bash
nanobot onboard
```

**3.配置** (`~/.nanobot/config.json`)

将以下两部分内容添加或合并到你的配置文件中（其他选项均使用默认值）。

*设置你的 API 密钥* ：
```json
{
  "providers": {
    "openrouter": {
      "apiKey": "sk-or-v1-xxx"
    }
  }
}
```

*设置你的模型*（可选择性指定服务方 —— 默认为自动检测）：
```json
{
  "agents": {
    "defaults": {
      "model": "anthropic/claude-opus-4-5",
      "provider": "openrouter"
    }
  }
}
```

**4. 聊天**

```bash
nanobot agent
```

## 消息通道

将纳米机器人连接到你常用的聊天平台。

| Channel | What you need |
|---------|---------------|
| **Feishu** | App ID + App Secret |

<details>
<summary><b>Feishu (飞书)</b></summary>

采用 **WebSocket** 长连接，无需公网 IP。

**1. 创建飞书机器人**
- 访问[飞书开放平台](https://open.feishu.cn/app)
- 创建新应用 → 启用**机器人**功能
- **权限**：添加 `im:message`（发送消息）和 `im:message.p2p_msg:readonly`（接收消息）
- **事件**：添加 `im.message.receive_v1`（接收消息）
- 选择**长连接模式**（需先运行 nanobot 建立连接）
- 在「凭证与基础信息」中获取**App ID**和**App Secret**


**2.配置**

```json
{
  "channels": {
    "feishu": {
      "enabled": true,
      "appId": "cli_xxx",
      "appSecret": "xxx",
      "encryptKey": "",
      "verificationToken": "",
      "allowFrom": ["ou_YOUR_OPEN_ID"]
    }
  }
}
```

> 长连接模式下，`encryptKey` 和 `verificationToken` 为可选配置。
> `allowFrom`：添加你的 open_id（向机器人发消息后可在 nanobot 日志中查看）。填写 ["*"] 表示允许所有用户。

**3. 运行**
```bash
nanobot gateway
```

</details>

<picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://api.star-history.com/svg?repos=zhiweihu1103/BreedingClaw&type=Date&theme=dark" />
    <source media="(prefers-color-scheme: light)" srcset="https://api.star-history.com/svg?repos=zhiweihu1103/BreedingClaw&type=Date" />
    <img alt="Star History Chart" src="https://api.star-history.com/svg?repos=zhiweihu1103/BreedingClaw&type=Date" />
</picture>
