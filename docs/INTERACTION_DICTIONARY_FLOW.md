# RXNRECer S3 交互字典流程文档

## 版本信息
- **版本**: 1.3.7
- **最后更新**: 2025-01-28
- **文档来源**: `data/dict/dict_rxnrecers3_prompt.json`

---

## 概述

本文档描述了RXNRECer S3模式中使用的交互字典（Prompt字典）及其工作流程。S3模式使用大语言模型（LLM）对S1和S2的预测结果进行推理、排序和解释。

---

## 交互字典结构

交互字典文件位置: `data/dict/dict_rxnrecers3_prompt.json`

字典包含4个不同的Prompt模板，用于不同的预测场景：

1. **prompt1**: 获取预测证据（Obtaining Evidence for Predictions）
2. **prompt2**: 细粒度结果排序和理由（Fine-grained Results Ranking and Justification for Proteins）
3. **prompt3**: 无UniProt ID的细粒度结果排序（Fine-grained Results Ranking and Justification for Proteins Without UniProt ID）
4. **prompt4**: 有UniProt ID的细粒度结果排序（Fine-grained Results Ranking and Justification for Proteins With UniProt ID）

---

## Prompt选择流程

### 流程逻辑

在 `rxnrecer/lib/llm/qa.py` 中的 `prompt_selector()` 函数负责选择使用的Prompt：

```python
def prompt_selector(rxnrecer_s1, rxnrecer_s2):
    set_s1 = set(rxnrecer_s1.split(cfg.SPLITER))
    set_s2 = set(rxnrecer_s2.split(cfg.SPLITER))
    
    if set_s1 == set_s2:
        return 'prompt1'  # S1和S2结果相同，使用prompt1
    else:
        return 'prompt2'  # S1和S2结果不同，使用prompt2
```

**选择规则**:
- 如果S1和S2的预测结果**完全相同** → 使用 `prompt1`
- 如果S1和S2的预测结果**不同** → 使用 `prompt2`

**注意**: 虽然字典中定义了 `prompt3` 和 `prompt4`，但在当前代码实现中主要使用 `prompt1` 和 `prompt2`。

---

## Prompt详细说明

### Prompt 1: 获取预测证据

**Prompt名称**: `Obtaining Evidence for Predictions`

**使用场景**: 当S1和S2预测结果相同时使用

**主要功能**:
- 为每个预测提供简洁的证据和解释
- 可参考UniProt、BRENDA、KEGG、ExPASy等权威来源
- 对非催化预测（"-"）提供直接证据

**关键考虑点**:

1. **注释来源**:
   - 如果有UniProt ID，从UniProt获取功能注释和EC编号作为基线证据
   - 如果没有UniProt ID，需要执行回退序列分析：
     - 模拟BLAST/DIAMOND同源性搜索，提出top 5同源物及其EC编号
     - 执行Pfam/InterProScan域检测，推断催化域和家族
     - 执行motif/活性位点扫描，识别关键保守残基或模式

2. **非催化预测**:
   - "-" 表示无催化活性
   - 提供直接证据（如缺少催化motif）

3. **多个反应**:
   - 如果预测了多个反应，为每个反应提供证据

4. **排序**:
   - 按可能性排序预测反应（1 = 最可能）

5. **置信度评分**:
   - 分配0.0-1.0的分数，反映序列特征对每个反应的支持强度
   - 如果选择"-"，给予最高置信度

6. **文献引用**:
   - 如果引用已发表的研究或数据库，引用真实、可验证的参考文献

**输出格式**:
```json
{
  "results": [
    {
      "reaction_id": "RHEA:XXXXX",
      "selected": "yes",
      "rank": 1,
      "confidence": 0.00-1.00,
      "reason": "简洁的证据摘要（无管道细节）"
    }
  ]
}
```

**规则**:
- 顶级对象必须只包含"results"键
- "results"必须是按可能性排序的对象数组（rank从1开始）
- 每个对象必须包含: `reaction_id`, `selected`, `rank`, `confidence`, `reason`
- 对于非催化预测，使用 `reaction_id: "-"`，`selected: "yes"`，并分配最高置信度
- 不要包含JSON之外的任何额外文本、解释或元数据

---

### Prompt 2: 细粒度结果排序和理由

**Prompt名称**: `Fine-grained Results Ranking and Justification for Proteins`

**使用场景**: 当S1和S2预测结果不同时使用

**主要功能**:
- 分析蛋白质序列和候选反应集合
- 确定哪些反应最可能由蛋白质催化
- 可参考UniProt、BRENDA、ExPASy等数据库

**关键考虑点**:

1. **注释来源**:
   - 如果有UniProt ID，从UniProt获取功能注释和EC编号
   - 如果没有UniProt ID，执行回退序列分析（同prompt1）

2. **非催化情况**:
   - "-" 表示预测缺乏催化活性
   - 如果所有反应都标记为"-"，或没有有效反应是合理的，选择"-"反应

3. **多个反应**:
   - 蛋白质可能有多个活性位点
   - 如果序列特征支持，选择多个反应

4. **排序**:
   - 为选定的反应分配基于可能性的排序（1 = 最可能）

5. **理由**:
   - 为每个决策提供简要解释
   - 引用保守域、活性位点残基、结构motif或已知催化机制

6. **置信度评分**:
   - 分配0-1之间的置信度分数
   - 基于蛋白质序列与已知催化特征的对齐程度
   - 如果没有有效反应是可能的，将最高置信度分配给"-"反应

**输出格式**:
```json
{
  "results": [
    {
      "reaction_id": "xxx",
      "selected": "yes" or "no",
      "rank": <integer>,     // 仅当 selected = "yes" 时
      "confidence": <float>,  // 0-1之间
      "reason": "基于序列和反应数据的解释"
    }
  ]
}
```

**规则**:
- 顶级对象必须只包含"results"键
- "results"必须是对象数组
- 对于 `selected: "yes"` 的条目，rank必须唯一且从1开始（最可能）
- 每个对象必须包含: `reaction_id`, `selected`, `confidence`, `reason`
- 仅当 `selected = "yes"` 时包含 `rank`
- 不要包含JSON之外的任何额外文本、解释或元数据

---

### Prompt 3: 无UniProt ID的细粒度结果排序

**Prompt名称**: `Fine-grained Results Ranking and Justification for Proteins Without UniProt ID`

**使用场景**: 当前代码中未直接使用，但可用于无UniProt ID的情况

**主要功能**: 与Prompt 2类似，但专门针对没有UniProt ID的蛋白质

**关键区别**:
- 自动执行BLAST/DIAMOND搜索
- 执行Pfam/InterProScan域扫描
- 执行motif扫描

---

### Prompt 4: 有UniProt ID的细粒度结果排序

**Prompt名称**: `Fine-grained Results Ranking and Justification for Proteins With UniProt ID`

**使用场景**: 当前代码中未直接使用，但可用于有UniProt ID的情况

**主要功能**: 与Prompt 2类似，但专门针对有UniProt ID的蛋白质

**关键区别**:
- 首先从UniProt检索功能注释和已知EC编号
- 使用UniProt数据作为基线证据

---

## S3工作流程

### 完整流程

```
输入: 蛋白质序列 + S1预测结果 + S2预测结果
  ↓
1. 选择Prompt (prompt_selector)
  ├─ S1 == S2 → prompt1
  └─ S1 != S2 → prompt2
  ↓
2. 构建查询字符串 (make_query_string)
  ├─ 提取UniProt ID和序列
  ├─ 获取反应详情 (RXN_details)
  └─ 构建JSON格式查询
  ↓
3. 批量处理 (make_query_batch)
  ├─ 读取Prompt字典
  ├─ 获取反应详情
  ├─ 为每行选择Prompt
  └─ 生成查询字符串
  ↓
4. LLM推理 (batch_chat / single_chat)
  ├─ 调用LLM API
  ├─ 传递系统Prompt和用户Prompt
  └─ 获取JSON格式响应
  ↓
5. 格式化输出 (format_rxn_output)
  ├─ 合并S3推理结果
  ├─ 格式化反应详情
  └─ 生成最终输出
```

### 代码实现位置

1. **Prompt选择**: `rxnrecer/lib/llm/qa.py::prompt_selector()`
2. **查询构建**: `rxnrecer/lib/llm/qa.py::make_query_string()`
3. **批量处理**: `rxnrecer/lib/llm/qa.py::make_query_batch()`
4. **LLM调用**: `rxnrecer/lib/llm/qa.py::batch_chat()` / `single_chat()`
5. **LLM客户端**: `rxnrecer/lib/llm/chat.py::Chat`
6. **输出格式化**: `rxnrecer/utils/format_utils.py::format_rxn_output()`

---

## 查询字符串构建

### 输入数据结构

```python
query_info = {
    "protein information": {
        "uniprot id": uniprot_id,
        "protein amino acid sequence": seq
    },
    "reaction information": {
        "predicted reaction 1": "RHEA:12345",
        "predicted reaction 2": "RHEA:67890",
        ...
    }
}
```

### 查询字典结构

```python
res_dict = {
    "PROMPT_SYS": system_prompt,  # 从字典中选择的prompt文本
    "PROMPT_USER": f"INPUT:\n{json.dumps(query_info, indent=2)}"
}
```

---

## LLM API调用

### Chat类接口

**类定义**: `rxnrecer/lib/llm/chat.py::Chat`

**初始化**:
```python
chat = Chat(
    name="openai/gpt-4.1",  # 模型名称
    url="https://api.openai.com/v1",  # API URL
    api_key="sk-...",  # API密钥
    proxy=None  # 可选代理
)
```

**调用方法**:
```python
response = chat.chat(
    message=dict_query['PROMPT_USER'],  # 用户消息
    system_prompt=dict_query['PROMPT_SYS'],  # 系统Prompt
    temperature=0,  # 温度参数
    max_tokens=10240,  # 最大token数
    response_format={"type": "json_object"}  # 响应格式
)
```

### 支持的LLM服务

- **OpenAI**: `https://api.openai.com/v1`
- **OpenRouter**: `https://openrouter.ai/api/v1`
- **Anthropic**: `https://api.anthropic.com`

### 环境变量配置

```bash
export LLM_API_KEY="your_api_key_here"
export LLM_API_URL="https://api.openai.com/v1"
```

---

## 输出格式

### S3模式输出结构

```json
{
  "reaction_id": "RHEA:24076",
  "prediction_confidence": 0.9999,
  "reaction_details": {
    "reactants": [...],
    "products": [...],
    "reaction_ec": "1.1.1.1",
    "reaction_equation": "...",
    "reaction_equation_ref_chebi": "...",
    "reaction_smiles": "..."
  },
  "reaction_rxnrecer_s3": {
    "selected": "yes",
    "ranking": 1,
    "confidence": 0.95,
    "reasoning": "基于序列特征和已知催化机制的解释..."
  }
}
```

### 字段说明

| 字段 | 类型 | 说明 |
|------|------|------|
| `reaction_id` | string | 反应ID (RHEA格式) |
| `prediction_confidence` | float | 预测置信度 (0-1) |
| `reaction_details` | object | 反应详细信息 |
| `reaction_rxnrecer_s3` | object | S3推理结果 |
| `selected` | string | 是否被选中 ("yes"/"no") |
| `ranking` | integer | 排序（1=最可能） |
| `confidence` | float | S3置信度 (0-1) |
| `reasoning` | string | 推理理由 |

---

## 错误处理

### 常见错误

1. **LLM API调用失败**:
   - 检查API密钥和URL配置
   - 检查网络连接
   - 查看错误日志

2. **JSON解析失败**:
   - LLM可能返回非标准JSON
   - 代码会捕获异常并返回默认结构

3. **Prompt选择错误**:
   - 确保S1和S2结果格式正确
   - 检查分隔符配置

### 错误处理代码

```python
try:
    response = chatcli.chat(...)
    content_str = response.choices[0].message.content
    records = json.loads(content_str)
    return records
except Exception as e:
    print(f"[ERROR] LLM query/parse failed: {e}")
    # 返回包含results字段的默认结构
    return {"results": [], "error": str(e)}
```

---

## 最佳实践

### 1. Prompt选择

- 确保S1和S2结果格式一致
- 使用标准分隔符（默认";"）

### 2. LLM配置

- 使用合适的模型（推荐GPT-4或Claude）
- 设置合理的max_tokens
- 使用JSON响应格式

### 3. 错误处理

- 始终检查LLM响应格式
- 提供默认值以防解析失败
- 记录错误信息用于调试

### 4. 性能优化

- 批量处理多个查询
- 使用缓存避免重复调用
- 合理设置批处理大小

---

## 配置参数

### 相关配置项

在 `rxnrecer/config/config.py` 中：

```python
# LLM API配置
LLM_API_URL = os.environ.get('LLM_API_URL', '')
LLM_API_KEY = os.environ.get('LLM_API_KEY', '')

# Prompt字典路径
FILE_DICT_RXNRECERS3_PROMPT = DATA_ROOT + 'dict/dict_rxnrecers3_prompt.json'

# 分隔符
SPLITER = ';'
```

---

## 扩展和定制

### 添加新Prompt

1. 在 `dict_rxnrecers3_prompt.json` 中添加新条目：
```json
{
  "prompt5": {
    "prompt_name": "Custom Prompt",
    "prompt_text": "..."
  }
}
```

2. 修改 `prompt_selector()` 函数以包含新逻辑

3. 更新 `make_query_string()` 以支持新Prompt

### 自定义LLM模型

修改 `batch_chat()` 中的默认模型：
```python
def batch_chat(..., llm_model='your-model-name'):
    ...
```

---

## 版本历史

- **v1.3.7**: 当前版本
  - 支持4种Prompt模板
  - 自动Prompt选择
  - JSON格式输出

---

## 参考资源

- **UniProt**: https://www.uniprot.org/
- **BRENDA**: https://www.brenda-enzymes.org/
- **KEGG**: https://www.genome.jp/kegg/
- **ExPASy**: https://www.expasy.org/
- **RHEA**: https://www.rhea-db.org/

---

## 联系信息

- **作者**: Zhenkun Shi
- **邮箱**: zhenkun.shi@tib.cas.cn
- **项目**: https://github.com/kingstdio/RXNRECer

