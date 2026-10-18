# PyMACS 技术推文大纲（定稿 v3）

> 三轮对抗性审查后定稿（含**论文全文 PDF 与补充材料 DOCX 的原文核对**）。所有技术断言均对照一手信源核实：
> 论文全文（*Eur. J. Med. Chem.* 2026, 316, 119038, DOI: 10.1016/j.ejmech.2026.119038 + 补充材料 mmc1.docx）、
> GitHub 仓库 `schurerlab/Pymacs`（README / docs / **实际脚本源码** / **MDP 模板** / GitHub API 文件清单）、官网 www.pymacs.com。
> 未核实的内容一律不写。
>
> 📌 使用说明：本文档分两区——**【正文结构】**即推文实际发布的内容（读者可见）；**【写作工作台】**为创作辅助材料（素材表 / 配图方案 / 修订说明），**不进入推文正文**。推文内部结构为：**第一部分 · 介绍（一~五章 + 参考资料）** + **第二部分 · Others（个人点评，放全文最后）**。

> 📷 封面头图（论文图形摘要）：`![](pyMACS_graph_abstract.jpg)` —— PyMACS 自动化流水线概念图

**拟定标题**：*PyMACS：把 GROMACS 从"手艺活"变成"流水线"——力场一致性、checkpoint 与 figurebook 的工程样本*

**备选标题**：*从 PDB 到 figurebook：PyMACS 如何用工程纪律驯服 GROMACS 的静默失败点*

**目标读者**：有 GROMACS 使用经验的计算化学/药物发现研究者 + 对"研究工具产品化"感兴趣的工程向读者

**核心叙事线索**：PyMACS 的价值不在"包装 gmx 命令"，而在三层工程决策——
① 力场家族一致性（CHARMM36 双力场 + CGenFF 配体 + 自动回退）；
② 运行状态可管理（checkpoint、状态文件、运行管理器）；
③ 分析闭环（轨迹 → 指标 → 网络 → figurebook 默认产出）。
并诚实交代边界：物理参数（MDP）仍是用户责任面，力场锁定限制适用社区。

**全文速览（3 句）**：
- 一条流水线：输入校验 → 拓扑 → 溶剂化 → 平衡 → 生产 MD → 轨迹分析 → figurebook，四步脚本闭环。
- 两种使用方式：新手交互问答 / HPC headless 脚本化，同一套参数两用。
- 三层工程脊梁：双力场自动回退、checkpoint 安全恢复、运行管理器一键提交（Triton/HPC）。

---

## 一、引言：MD 工作流的"暗面"

### ① 被低估的工程问题

- MD 领域的叙事重心在物理模型与采样，但**工作流工程**长期被视为"已解决"
- 现实：一个标准蛋白-配体体系的 GROMACS 设置涉及 10+ 手工步骤，每一步都有静默失败的入口
- 论文原话（摘要）：*"errors at any stage can compromise results and discourage non-expert users"*
- **硬素材**：PyMACS 官方 troubleshooting 列出 6 类 grompp 失败原因——力场不匹配、重复 include、CGenFF 版本漂移、残基名不匹配、缺力场目录、缺 checkpoint 文件。工程问题之具体，文档自己就是证据

> 📷 配图（论文 Fig 1 · 典型 MD 工作流）：`![](pyMACS_f1.jpg)` —— 准备/参数化 → 模拟 → 分析/可视化，直观呈现"手工步骤多"

### ② 生态谱系：PyMACS 的定位

> ⚠️ 写法要求：本节为**作者整理的生态对照**，各工具描述挂独立信源（引用编号兜底）；以客观陈述呈现，不写成"论文中的谱系"。

- **gmxapi**（GROMACS 官方 Python 接口，Irrgang et al., PLoS Comput Biol 2022）：底层 API，不面向完整工作流
- **GromacsWrapper / gromacs_py**：thin wrapper 族，简化 gmx 调用、可建拓扑，但无配体参数化闭环、无分析默认产出
- **StreaMD**（Ivanova et al., J Cheminform 2024）：高通量端到端流水线，面向批量化而非精细交互控制
- **LLM-agent 族**：GROMACS-Agent（ChemRxiv 预印本）、MDCrow（2026，~29 引）、gromacs_copilot——把"工作流生成"交给 LLM，可靠性边界尚未建立
- **PyMACS 的差异点**：在 CHARMM36/CGenFF 家族一致性承诺下，提供"输入校验 → 拓扑 → 模拟 → 分析 → figurebook"闭环，同时把 MDP、box、cutoff 的控制权留给用户

> **推文锚点**：谱系对比表，横轴"分析闭环深度"、纵轴"力场一致性承诺"，标出 PyMACS 位置。

---

## 二、架构：三层设计，一图看懂

*本节论点：PyMACS 不是单一 CLI，而是"脚本链 + 双模式 + 运行管理器"的分层架构，每层解决一类问题。*

### 2.1 三层结构总览

- **脚本链层**：1（setup）→ 2（EM/NVT/NPT/生产）→ 3A / 3B / 3C / 3_PROTAC（分析）→ 4（figurebook）；每步可独立复用（如只用 3A 分析已有轨迹）；另有辅助脚本 `2A_AutoGMXrestart.py`（生产恢复+后处理）、`3A_..._chunks`（分块分析变体，默认每块 1000 帧）
- **双模式层**：interactive Q&A（新手引导）vs `--headless` + 全 flag 集（HPC / tmux / 脚本化）；同一套参数两用，避免"教程一套、生产一套"的割裂
- **运行管理器层**：`pymacs_run.py` 提供 7 个子命令——`configure`（问答写 `pymacs_run.json`）→ `validate` → `setup`（驱动 MPI 版 Step 1）→ `submit`（渲染并提交 LSF 作业，可选 `--with-analysis` 链式提交分析作业）→ `resume`（checkpoint 续跑作业）→ `submit-analysis` → `status`（查看运行状态）；HPC 档案在 `hpc_profiles/`（`triton.json` + 两个 LSF 模板 + `leela/`）；所有执行命令写入 `pymacs_commands.log`
- **技术讨论点**：为什么"脚本链 + 管理器"而非单一 CLI？——降低模块耦合、允许中途介入；每步产出即中间品，契合 run-folder 工作纪律

> 📷 配图（论文 Fig 2 · Script 1 系统准备与参数化）：`![](pyMACS_f2.jpg)` —— 脚本链层的"第一步"长什么样

### 2.2 可靠性工程

- **Checkpoint 安全恢复**：`md_0_1.tpr / md_0_1.cpt` 为核心；`--resume --production_only` 跳过 EM/NVT/NPT 只续生产；`--force_restart` 强制重跑；细节坑：restart 时 GROMACS 不读 md.mdp，`--ns` 必须 ≤ 原 tpr 长度
- **恢复 + 后处理自动化（2A 脚本）**：自动续跑（有 .cpt 则续、无则新跑）→ GPU 失败自动回退（检测 virtual site/update 不兼容时切 CPU-update 混合模式）→ 轨迹 wrap/center（Protein+Ligand）→ 提取结合口袋（5.0 Å，每 5 帧采样）
- **状态文件与可复现**：链命名缓存在 **`.chainmap.json`**（论文 §2.1 / Table 2 明确记载，用于跨运行一致命名）；组分溯源用 `atomIndex.txt` + `pymacs_system_registry.json` / `pymacs_components.json`；`pymacs_run.json`（运行决策永久记录，可 review/编辑）、`pymacs_commands.log`（命令审计）、`pymacs_analysis_cache.py`（分析缓存）；论文 §3.1：力场、末端、链名、盒子、MDP 来源、索引、执行命令**全部记录在日志**；运行纪律：overwrite-by-default + `RUNS/System_run01` 目录模式
- **资源工程**：`--compute auto` 自动 GPU 探测 → 多 GPU 交互询问 → 无 GPU 时 CPU 回退；EM 阶段刻意避开 PME-on-GPU（非动力学最小化器不兼容）；MPI 变体族（`PYMACS_GMX_BIN`、`--external-mpi --mpi-ranks --mpi-launcher`、srun/mpirun、SLURM/LSF）

> 📷 配图（论文 Fig 3 · Script 2 平衡与生产 MD）：`![](pyMACS_f3.jpg)` —— EM → NVT → NPT → 生产，checkpoint 落点

---

## 三、核心实现：力场、配体与分析闭环

*本节论点：PyMACS 的"窄而深"体现在力场家族一致性；它的"广"体现在分析产出全流程覆盖。*

### 3.1 力场家族与配体参数化

- **双力场 + 自动回退**：内置 `charmm36.ff`（经典）与 `charmm36_ljpme-jul2022.ff`（LJ-PME 变体）；部分体系（尤其 CGenFF 配体体系）会先试其一、失败自动换另一个——容错优先的务实设计
- **论文的力场上下文（§1.2，写作时带出处）**：标准 CHARMM36 工作流常用 ~10–12 Å force-switching 处理 Lennard-Jones 相互作用；CHARMM36/LJ-PME 对长程 LJ 采用 particle-mesh-Ewald 处理，可用 9 Å 非键 cutoff 提升效率——注意这是论文对力场背景的陈述，PyMACS 的实际默认值仍由用户可替换的 MDP 模板决定
- **论文默认体系**：SPC 水模型、0.15 M NaCl 中和、**默认立方盒 + 1.0 nm buffer**（向后兼容，论文 §2.1；示例常改用 dodecahedron 省溶剂）
- **诚实边界（必写）**：内置 MDP 模板为教程级通用值（EM/NPT/MD：rvdw/rcoulomb 1.0 nm；NVT：1.5 nm；DispCorr；无 force-switch），即"可用但非 CHARMM 优化默认"；论文明确"默认 MD 参数定义在外部 .mdp 模板、不硬编码进脚本"（§2.1）；精确协议靠 `--mdp-dir / --em-mdp / --nvt-mdp / --npt-mdp / --md-mdp` 用户替换 + `--equilibration-plan` JSON 多级平衡
- **CGenFF 配体参数化（四路径）**：① 导入现成 GROMACS 文件（.itp/.prm，BYOP）② CGenFF server 输出 `.str` + `.cgenff.mol2` → `cgenff_charmm2gmx_py3_nx2.py` → GROMACS include（推荐）③ 本地 SILCSBio 全自动（提取→生成→转换→整合）④ **远程 CGenFF funnel 模式**（token 鉴权的 HTTPS 服务，适合无本地 CGenFF 的 HPC 环境）
- **版本锁定**：验证于 CGenFF 4.6；新版输出可能不兼容——文档明示的兼容边界
- **共价配体与核酸**：LINK/CONECT 检出时**默认警告或停止**（论文 §2.1 原文），`--allow-covalent-ligand` 专家覆盖；核酸残基**从配体自动检测中排除**（论文 §2.1），走 `--mode biological` + `--nucleic-acid-termini` / `--chain-types`

### 3.2 分析引擎与 figurebook

- **基础指标**：RMSD、RMSF、Rg、DSSP 二级结构演化（SSE fractions、helix/sheet persistence）

> 📷 配图（论文 Fig 4 · Script 3A/3B/4 轨迹分析与可视化）：`![](pyMACS_f4.jpg)` —— 分析闭环的骨架图
- **配体中心**：pocket（默认 5.0 Å）、contact（默认 4.0 Å）、氢键分类（3.5 Å / 135°）、`binding_importance_ranking` 热点残基排名、`--min_contact_frac`（默认 0.10）控制网络拓扑；论文明示疏水/离子分类属"筛选级注释，非自由能替代"（§2.3）
- **Interface RIN（3C）**：无配体多链体系自动分发；chain-pair 接触持久性、residue-pair 热图、链级网络；`--interface-chain-pairs A:B`、`--interface-frame-step`、`--interface-max-edges`
- **PROTAC（3_PROTAC）**：组件级 RMSD/RMSF、三元接触、几何、水桥（3.5 / 6.0 Å）、QC manifest、`--quick-test` 冒烟测试
- **跨运行对比**：`00_Pymacs_X_analysis.py`（可选，Plotly HTML 仪表盘）——论文 §2.4 确认：扫描完成运行目录，生成 apo vs ligand / ligand vs ligand / 受体家族 / 跨物种静态 HTML 仪表盘 + CSV
- **Figurebook（4）**：`4_MDfigs.txt` 控制收录与顺序，`4_GraphNotes.txt` 加注；PROTAC 模式自动切换 manifest（`protac_figure_manifest.csv`）
- **技术讨论点**：figurebook 是"分析闭环"的收口——对比传统"跑 gmx → 导 xvg → matplotlib → 排版"，PyMACS 把它变成默认产出；但注意 figurebook 不产生新科学结论，只做叙事打包

> 📷 配图（论文 Fig 5 · 代表性分析输出六面板）：`![](pyMACS_f5.jpg)` —— (A) 蛋白 RMSD (B) 配体 RMSD (C) Rg (D) 接触图 (E) 热点排名 (F) NETWORX 网络，推文里最直观的"产出长这样"

---

## 四、实战：最小可运行示例

*本节论点：让读者 30 分钟内能亲手跑通一次，并知道生产级配置的旋钮在哪。*

### 4.1 四步跑通最小示例

- 前置：GROMACS 2022+（companion installer：Joey305/gromacs-installation）+ 两个 conda 环境（`environment_cgenff.yml` / `environment_mdanalysis.yml`，激活名 cgenff / mdanalysis）
- 输入：PDB / CIF / mmCIF，或 `--fetch-pdb 9UWJ`；`--box-type dodecahedron --box-distance 1.0`
- 命令块（README 快速上手提炼）：

```bash
# Step 1 setup（结构清理、拓扑、溶剂化、离子、EM 输入）
conda activate cgenff
python 1_AutomateGromacs.py --pdb complex.pdb --ligand A1D

# Step 2 平衡 + 生产（--ns 0.25 快速测试；--compute auto 自动 GPU）
conda activate mdanalysis
python 2_AutomateGromacs.py --ns 0.25

# Step 3 分析（配体模式自动产出 RMSD/RMSF/contacts）
python 3A_AutomateGromacs.py

# Step 4 figurebook PDF
python 4PDF4MD.py
```

- 仓库内置 4 个可复现示例：Example1（CPD32_9G94 + A1D，蛋白-配体）、Example2（9UWJ/AVPR1A + A1E + CLR，配体/辅因子）、Example3（1URN，RNA/蛋白生物装配）、Example4（5T35，PROTAC 三元复合物）；大文件走 Git LFS
- **实测数据（补充材料 S2/S8，可当"真实跑分"素材）**：A1D 示例——RTX 3060 + 12 线程，EM 364 步收敛、NVT/NPT 各 5 万步、生产 250 万步（5 ns），354.959 ns/day；1URN RNA 示例 0.25 ns，106.954 ns/day

> 📷 配图（补充材料 SI，已保存为 si_s1/si_s2.png）：`![](si_s1.png)` 交互式 setup 运行截图；`![](si_s2.png)` 交互式 simulation 运行截图（带 354.959 ns/day 真实跑分）——"这不是演示，是真跑出来的"

### 4.2 生产配置、恢复与批量

- **MDP 与协议控制**：`--mdp-dir` 系列替换模板、`--index-file`、`--equilibration-plan` JSON 多级平衡、`--ns` 生产时长、restraints
- **中断恢复**：checkpoint 细节（见 §2.2）+ `--ns` 的坑；HPC 上 `pymacs_run.py resume --profile triton` 一键提交续跑作业
- **批量与对比**：3A 按单运行目录分析（"batch-friendly"为设计目标）；多运行对比走 `00_Pymacs_X_analysis.py`（论文 §2.4：apo vs ligand / ligand vs ligand / 受体家族 / 跨物种）

---

## 五、评估：优势、边界与启示

*本节论点：评价一款科研工具，要看它"替用户记住了什么"，也要看它"把什么留给了用户"。*

### 5.1 优势

- **力场家族一致性**：从拓扑生成到分析 cutoff 都在 CHARMM 框架内，减少"参数混用"
- **分析闭环 + 产品化工程**：GPU 自适应、checkpoint、运行管理器、HPC 一键提交、命令审计日志、测试目录（tests/ 含 run 管理器/缓存/funnel 测试）
- **可复现接口**：registry 文件、pymacs_run.json、MDP 模板外置——"可复现"是操作化的，不是口号

### 5.2 边界

- **力场锁定** CHARMM 家族：AMBER/OPLS/GROMOS 用户需自带 MDP 与转换
- **MDP 内置模板为教程级**，生产协议仍需用户判断——诚实的边界，不是 bug
- **依赖 GROMACS 预装**：companion installer 缓解，但编译与 GPU 配置仍是门槛
- **Docking / pose 生成不在范围内**：论文 §2.1 明示"假定结合 pose 已由上游生成，PyMACS 自身不做 docking 或 pose 预测"；消费上游结果（AlphaFold / RoseTTAFold / Boltz-2 / GNINA / Chai-1 / DiffDock / AutoDock Vina 等）
- **方法论范围（论文 §3.3 原文）**：不做自由能计算、增强采样、电子极化；不做 pose 枚举/排序、质子化/电离/互变异构枚举；RNA 支持以 setup/模拟为主，RNA-ligand 基准与 RNA 特异性分析是未来工作
- **CGenFF penalty score**：论文 §1.2 定位为"帮助判断是否需要额外参数优化"的信号，无自动专家反馈机制，官网明示高 penalty 需人工精修
- **LJ-PME 参数敏感性**：论文 §3.3 提示脂质体系需专门优化、蛋白-配体体系需按靶点单独验证
- **LLM-agent 族对比**：规则流水线 vs 生成式工作流，可靠性与灵活性的边界不同

### 5.3 启示

- **"替代的是流程记忆，不是物理判断"**：PyMACS 把资深用户的流程内隐知识（步骤顺序、恢复策略、分析默认阈值）编码进默认值，但把物理参数旋钮留给用户——自动化是在正确的地方暴露正确的旋钮
- 可复制性：这种"力场家族 + 分析闭环 + HPC 提交"模式能否在 OpenMM / AMBER 生态复制？

**参考资料（归入第一部分末尾，GB/T 7714 顺序编码 [1]–[8]；正文按首次出现编号、条目含 URL）**

---

## 第二部分 · Others：个人点评（全文最后，由作者本人撰写）

> 第一部分（介绍：一~五章 + 参考资料）到此为止。**Others 是全文最后一部分**（放在参考资料之后），由作者本人撰写个人点评；**结尾祝福也写在这里**。以下提供干条条素材，可自由改写、删减或增补。

【干条条素材】
① 最值钱的不是算法，是工程纪律——无新力场/新采样/新模型，价值在把流程做扎实。
② 方法学上没有、工程上有——核心概念"工作流产品化"；押注"确定性自动化"，可当 LLM 生成式路线的对照基线；金句：把领域最佳实践做成默认值，把可复现性做成机制，把流水线做成产品。
③ 为什么能发 EJMC——"Comp small drug design"特刊定位 + 受众精准（药化/结构生物）+ 资源完整（代码/示例/补充材料/开放获取/50 篇参考文献）。
④ 保留意见（四条）——MDP 教程级、力场锁定、RNA 分析偏弱、依赖 GROMACS 预装；别当"自动写好 MD"的黑盒。
⑤ 适用人群——CHARMM36/CGenFF 用户当标准模板评估；缺算力人力的药化小组最受益；AMBER 生态借鉴架构；LLM 爱好者当对照基线。
⑥ 收束金句（可选）——"把流程做对，本身就是创新；它替我们记住步骤，把判断留给我们。"

【结尾祝福】由作者写在此处（技能风格示例：祝各位好~ / 科研顺利~ / 欢迎讨论和交流）。

---

# 写作工作台（创作辅助材料，不进入推文正文）

> 以下全部为创作辅助材料（素材表 / 配图方案 / 修订说明），推文读者不会看到这些。

**📌 写作总则（作者定，适用于今后所有推文）**：
- 推文分两大类：**① 文章/工具介绍**；**② 原创教程**（针对某一具体事情）。
- 两类统一采用**两段式结构**：
  - **第一部分 · 客观陈述**（助手撰写）：第三方视角、客观、深入、有条理；不掺"笔者怎样怎样"等主观表述；参考文献归本部分。
  - **第二部分 · Others**（作者本人撰写）：个人点评/感受 + 结尾祝福，作者自由发挥；助手只提供干条条素材。
- 助理的工作范围 = 第一部分 + 提供第二部分素材；第二部分正文由作者完成。

## 附录 A：推文配套素材表

| 素材 | 用途 | 来源 | 状态 |
|---|---|---|---|
| 三层架构图（脚本链 / 双模式 / 运行管理器 7 子命令） | 第二节架构说明 | 基于 README + pymacs_run.py 源码重绘 | 需自制 |
| 力场策略说明（双文件夹 + auto-fallback + MDP 边界） | 第三节力场讨论 | README + 仓库 MDP 模板 | 已核实 |
| CGenFF 四路径流程图（含远程 funnel） | 第三节配体参数化 | README + pymacs_run.py / TRITON_CGENFF_FUNNEL.md | 已核实 |
| Example1 figurebook PDF | 第三节分析引擎示例 | GitHub `Example_Choices/Example1/completed_run/`（Git LFS） | 可直接引用 |
| 最小四步命令块 | 第四节实战 | README quick start 提炼 | 已核实 |
| 实测运行数据（S2/S8：354.959 / 106.954 ns/day 等） | 第四节实战彩蛋 | 补充材料 mmc1.docx | 已核实 |
| 生态谱系对比表 | 第一节/第五节 | 需按独立信源自制（gmxapi / GromacsWrapper / gromacs_py / StreaMD / MDCrow…） | 信源清单已备 |

## 附录 B：推文配图方案

> 正文已插入图片引用（md 语法）；图与大纲同在 `PyMACS/` 目录，相对路径直接可用。SI 的 S1/S2 已从 mmc1.docx 提取为 `si_s1.png` / `si_s2.png`；其余按需提取（见下）。

| 推文位置 | 图 | md 引用 | 图例要点 | 来源 | 状态 |
|---|---|---|---|---|---|
| 头图/封面 | 图形摘要 | `![](pyMACS_graph_abstract.jpg)` | PyMACS 自动化流水线概念图 | 本地 | ✅ 已就位 |
| §1 引言① 工程问题 | 论文 Fig 1 | `![](pyMACS_f1.jpg)` | 典型 MD 工作流：准备 → 模拟 → 分析/可视化 | 本地 | ✅ 已就位 |
| §2.1 脚本链层 | 论文 Fig 2 | `![](pyMACS_f2.jpg)` | Script 1 系统准备与参数化流程 | 本地 | ✅ 已就位 |
| §2.2 可靠性工程 | 论文 Fig 3 | `![](pyMACS_f3.jpg)` | Script 2：EM → NVT → NPT → 生产 MD | 本地 | ✅ 已就位 |
| §3.2 分析闭环 | 论文 Fig 4 | `![](pyMACS_f4.jpg)` | Script 3A/3B/4：轨迹分析 → 网络 → figurebook | 本地 | ✅ 已就位 |
| §3.2 分析产出 | 论文 Fig 5 | `![](pyMACS_f5.jpg)` | 代表性输出六面板（RMSD/Rg/接触图/热点排名/NETWORX） | 本地 | ✅ 已就位 |
| §4.1 实战·setup | SI S1 | `![](si_s1.png)` | 交互式 setup 运行截图（CPD32_9G94 + A1D） | 已从 docx 提取 | ✅ 已保存 |
| §4.1 实战·模拟 | SI S2 | `![](si_s2.png)` | 交互式 simulation 截图（RTX 3060，354.959 ns/day） | 已从 docx 提取 | ✅ 已保存 |
| §3.2 figurebook | SI S4 | `![](si_s4.png)` | figurebook 生成运行截图 | docx 内嵌 image4 | ⏳ 待保存 |
| §5.2 多体系能力（可选） | SI S7 或 S9 | `![](si_s7.png)` | RNA 1URN setup 截图（或 PROTAC S9） | docx 内嵌 image7/9 | ⏳ 待保存（可选） |

**备用未列入**：S3（分析截图）、S5/S6（cofactor 体系）、S8（RNA 模拟）、S10（PROTAC 模拟）——写作时按篇幅取舍。

**SI 图提取**：mmc1.docx 本质是 zip，10 张图位于其 `word/media/image1.png … image10.png`（顺序对应 S1–S10）。已完成：image1→`si_s1.png`、image2→`si_s2.png`。如需 S4/S7（figurebook 示例、RNA 体系），按同法提取命名后放入本目录即可。

---

## 附录 C：修订说明 v3（三轮对抗性审查）

**第一轮（事实核查）**：
1. ~~删除 `.chainmap.json`~~——**已被第三轮推翻**：论文 §2.1/Table 2 明确记载，已恢复（与 README 的 atomIndex.txt/registries 并存）。
2. 删除未核实数值（LJ-PME 9 Å / 10–12 Å force-switching）——**第三轮确认论文 §1.2 确有表述**（属力场上下文陈述，非 PyMACS 默认值）；**已按用户确认写入正文 §3.1 并标注论文出处，同时注明"非默认值、默认由 MDP 模板决定"**。
3. "力场承诺"降级为"家族一致性 + 边界诚实"——内置 MDP 为教程级通用模板。
4. 谱系表重定位为作者独立整理，补 MDCrow、gmxapi。
5. 分析引擎扩写：DSSP、binding-importance ranking、水桥、Interface RIN、跨运行对比。

**第二轮（GitHub API + 源码级核实）**：
6. ✅ 确认 `pymacs_run.py`（27 KB）、`hpc_profiles/`（triton.json + 2 个 LSF 模板 + leela/）、Triton/Leela 环境 YAML、`00_Pymacs_X_analysis.py`、`tests/` 全部真实存在于 main 分支——运行管理器卖点成立；子命令实为 **7 个**（含 submit-analysis / status），并驱动 **chunked 分析**（3A_chunks_MPI，--chunk-size 1000）。
7. 🆕 新增发现并纳入：`2A_AutoGMXrestart.py`（恢复+后处理，GPU 失败回退、口袋提取）、CGenFF **远程 funnel 第四路径**、`pymacs_analysis_cache.py`、`pymacs_commands.log` 命令审计、4 个示例 PDB 与目录。
8. 措辞修正：备选标题去掉"10 个"具体数字。（共价配体措辞见第 9 条，第三轮按论文原文改回）

**第三轮（论文全文核对：pypdf 提取 PDF 正文 + python-docx 读补充材料 mmc1.docx）**：
9. ✅ `.chainmap.json` 恢复（见第 1 条）；共价配体措辞改回论文原文"默认警告或停止"。
10. ✅ 补论文默认体系：SPC 水、0.15 M NaCl、默认立方盒 1.0 nm buffer；补 §3.3 方法论范围（free-energy / 增强采样 / 极化 / pose 枚举均不在范围内）、RNA 分析边界、LJ-PME 脂质敏感性；§3.1"默认 MD 参数外置于 .mdp、不硬编码"获论文背书。
11. ✅ 对比仪表盘获论文 §2.4 背书；补充材料 S1–S10 实测数据（EM 364 步、354.959 / 106.954 ns/day 等）入素材表。
12. 可选彩蛋：论文披露使用 Deep Research（ChatGPT）辅助大纲与配图（原写入"结语"；架构调整后该信息并入 Others 素材 ③）。
13. 新增"六、我们的看法"观点板块：亮点（工程纪律而非算法）/ 创新定位（工作流产品化，非范式革命）/ 选刊逻辑（EJMC 特刊 + 受众精准）/ 保留意见 / 适用人群——作者观点，标注非论文立场。（后被 item 15 移除，观点转为 Others 干条条素材）
14. 按用户意见再调整：② 改为"方法学上没有、工程上有"（不谈范式，核心概念改"工作流产品化"+"确定性自动化"）；原结语章并入 ⑥ 收束；文档分两区（正文结构 / 写作工作台），附录明确标注不进入推文。
15. ✅ 最终架构（用户拍板）：全文分两部分——**第一部分 · 介绍（一~五章 + 参考资料，参考资料归介绍）**；**第二部分 · Others（个人点评，放全文最后，位于参考资料之后）**，Others 内容与结尾祝福由作者本人撰写，助手仅提供干条条素材；"六、我们的看法"从正文移除，改为 Others 素材；大纲文档分区命名改为【正文结构】/【写作工作台】以避开与推文内部"第二部分"混淆。
16. ✅ 第一部分全面客观化（用户要求）：删除全部"笔者"主观表述（开头自我介绍、"笔者很认同/觉得/愿意"等）、感叹与 emoji；改为第三方客观文风（"论文明确指出""一个值得注意的启示是"）；此安排写入"写作总则"，适用于今后所有推文。

> 数值出处约定：标注"论文"者出自正文 §1.2 / §2.1 / §2.3 / §2.4 / §3.3 / Table 2 / 补充材料 mmc1.docx；其余出自仓库 README / MDP 模板 / 脚本源码（GitHub API 与 raw 核实）。
