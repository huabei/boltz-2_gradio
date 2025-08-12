import gradio as gr
import yaml
import subprocess
import tempfile
import shutil
import os
import json
import datetime
from pathlib import Path
import base64 # 新增导入

def get_available_gpus():
    """检测可用的GPU数量"""
    try:
        result = subprocess.run(['nvidia-smi', '--query-gpu=index', '--format=csv,noheader,nounits'], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            gpu_lines = [line.strip() for line in result.stdout.strip().split('\n') if line.strip()]
            gpu_count = len(gpu_lines)
            return gpu_count if gpu_count > 0 else 1
        else:
            return 1
    except (subprocess.TimeoutExpired, FileNotFoundError, Exception):
        return 1

def create_formatted_affinity_markdown(affinity_data):
    """将亲和力JSON数据格式化为易于阅读的Markdown。"""
    if not affinity_data:
        return "未生成亲和力数据。"

    value = affinity_data.get("affinity_pred_value", "N/A")
    prob = affinity_data.get("affinity_probability_binary", "N/A")

    # 格式化数值
    if isinstance(value, float):
        value_str = f"{value:.2f}"
    else:
        value_str = "N/A"
        
    if isinstance(prob, float):
        prob_str = f"{prob:.2%}"
    else:
        prob_str = "N/A"

    # 添加来自文档的解释
    explanation = (
        "**解读亲和力分数:**\n"
        "- **预测亲和力值 (Predicted Affinity Value)**: 以 `log(IC50)` 形式报告，其中 IC50 单位为 `μM`。**值越低，预测的结合能力越强。**\n"
        "  - 例如: -3 (强结合剂), 0 (中等结合剂), 2 (弱结合剂/诱饵)。\n"
        "- **结合概率 (Binding Probability)**: 范围从 0 到 1，表示配体是结合剂的预测概率。"
    )
    
    md = f"""
### 亲和力预测结果
| 指标 | 预测值 |
| :--- | :--- |
| **预测亲和力值** | **`{value_str}`** |
| **预测亲和力(Kcal/mol)** | **`{((6-value) * 1.364):.2f}`** |
| **结合概率** | **`{prob_str}`** |

---
{explanation}
    """
    return md

def create_formatted_confidence_markdown(confidence_data):
    """将置信度JSON数据格式化为易于阅读的Markdown。"""
    if not confidence_data:
        return "未生成置信度数据。"

    score = confidence_data.get("confidence_score", "N/A")
    iptm = confidence_data.get("iptm", "N/A")
    plddt = confidence_data.get("complex_plddt", "N/A")

    if isinstance(score, float): score_str = f"{score:.3f}" 
    else: score_str = "N/A"
    if isinstance(iptm, float): iptm_str = f"{iptm:.3f}"
    else: iptm_str = "N/A"
    if isinstance(plddt, float): plddt_str = f"{plddt:.3f}"
    else: plddt_str = "N/A"

    explanation = (
        "**解读置信度分数 (范围 0-1, 越高越好):**\n"
        "- **综合置信度 (Confidence Score)**: 用于对模型排序的聚合分数 (0.8 * complex_plddt + 0.2 * iptm)。\n"
        "- **iptm**: 预测的界面TM-score，衡量链间相互作用预测的准确性。\n"
        "- **complex_plddt**: 复合物的平均pLDDT，衡量原子位置预测的局部置信度。"
    )

    md = f"""
### 结构置信度分数
| 指标 | 预测值 |
| :--- | :--- |
| **综合置信度** | **`{score_str}`** |
| **iptm (链间)** | **`{iptm_str}`** |
| **complex_plddt (局部)** | **`{plddt_str}`** |

---
{explanation}
    """
    return md
# 添加新的 get_molstar_html 函数
def get_molstar_html(mmcif_base64):
    return f"""
    <iframe
        id="molstar_frame"
        style="width: 100%; height: 600px; border: none;"
        srcdoc='
            <!DOCTYPE html>
            <html>
                <head>
                    <script src="https://cdn.jsdelivr.net/npm/@rcsb/rcsb-molstar/build/dist/viewer/rcsb-molstar.js"></script>
                    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/@rcsb/rcsb-molstar/build/dist/viewer/rcsb-molstar.css">
                </head>
                <body>
                    <div id="protein-viewer" style="width: 1200px; height: 400px; position: center"></div>
                    <script>
                        console.log("Initializing viewer...");
                        (async function() {{
                            // Create plugin instance
                            const viewer = new rcsbMolstar.Viewer("protein-viewer");

                            // CIF data in base64
                            const mmcifData = "{mmcif_base64}";

                            // Convert base64 to blob
                            const blob = new Blob(
                                [atob(mmcifData)],
                                {{ type: "text/plain" }}
                            );

                            // Create object URL
                            const url = URL.createObjectURL(blob);

                            try {{
                                // Load structure
                                await viewer.loadStructureFromUrl(url, "mmcif");
                            }} catch (error) {{
                                console.error("Error loading structure:", error);
                            }}
                      }})();
                    </script>
                </body>
            </html>
        '>
    </iframe>"""

def get_initial_molstar_html():
    """为 Gradio 界面生成初始的 Mol* 查看器 HTML。"""
    default_cif_path = Path("example.cif") # 假设 example.cif 在应用根目录
    if default_cif_path.exists():
        try:
            mmcif_bytes = default_cif_path.read_bytes()
            mmcif_base64 = base64.b64encode(mmcif_bytes).decode('utf-8')
            return get_molstar_html(mmcif_base64)
        except Exception as e:
            error_message = f"加载默认 example.cif 时出错: {e}"
            print(error_message)
            return f"<div style='height: 600px; display: flex; align-items: center; justify-content: center;'><p>{error_message}</p></div>"
    else:
        # 如果 example.cif 不存在，返回一个空的 Mol* 查看器或提示信息
        return get_molstar_html("") # 传递空字符串，让 Mol* 内部处理或显示无数据

def run_boltz_prediction(
    sequences_config,
    use_msa_server,
    use_potentials,
    recycling_steps, 
    diffusion_samples,
    enable_affinity_prediction,
    affinity_binder_id,
    gpu_count
):
    """
    一个完整的函数，用于生成YAML，运行Boltz，并处理输出。
    支持多种分子类型和多条链预测。
    """
    # 1. 输入验证
    if not sequences_config or len(sequences_config) == 0:
        return "错误：至少需要添加一个分子序列。", None, None, None, None, None, None
    
    # 验证亲和力预测设置
    if enable_affinity_prediction:
        if not affinity_binder_id.strip():
            return "错误：启用亲和力预测时必须指定结合分子的链ID。", None, None, None, None, None, None
        
        # 检查结合分子ID是否存在于配置中
        chain_ids = [seq.get("chain_id", "").strip() for seq in sequences_config if seq.get("chain_id", "").strip()]
        if affinity_binder_id.strip() not in chain_ids:
            return f"错误：结合分子链ID '{affinity_binder_id}' 不存在于当前分子列表中。", None, None, None, None, None, None

    # 创建tmp目录来存放所有文件（便于查看和调试）
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    tmp_base_dir = Path("tmp")
    tmp_base_dir.mkdir(exist_ok=True)
    run_dir = tmp_base_dir / f"boltz_run_{timestamp}"
    run_dir.mkdir(exist_ok=True)
    
    initial_3d_html = "<div style='height: 600px; display: flex; align-items: center; justify-content: center;'><p>等待预测结束...</p></div>"

    try:
        # 2. 生成 YAML 配置文件
        input_dir = Path(run_dir) / "input"
        output_dir = Path(run_dir) / "output"
        input_dir.mkdir()
        output_dir.mkdir()

        config_name = "prediction_config"
        yaml_path = input_dir / f"{config_name}.yaml"
        
        # 构建YAML序列部分
        yaml_sequences = []
        for seq_config in sequences_config:
            chain_id = seq_config["chain_id"].strip()
            mol_type = seq_config["mol_type"]
            sequence = seq_config["sequence"].strip()
            
            if not chain_id or not sequence:
                continue  # 跳过空的配置
            
            if mol_type == "蛋白质":
                yaml_sequences.append({
                    "protein": {
                        "id": chain_id,
                        "sequence": sequence
                    }
                })
            elif mol_type == "DNA":
                yaml_sequences.append({
                    "dna": {
                        "id": chain_id,
                        "sequence": sequence
                    }
                })
            elif mol_type == "RNA":
                yaml_sequences.append({
                    "rna": {
                        "id": chain_id,
                        "sequence": sequence
                    }
                })
            elif mol_type == "配体(SMILES)":
                yaml_sequences.append({
                    "ligand": {
                        "id": chain_id,
                        "smiles": sequence
                    }
                })
            elif mol_type == "配体(CCD)":
                yaml_sequences.append({
                    "ligand": {
                        "id": chain_id,
                        "ccd": sequence
                    }
                })

        # 构建完整的YAML结构
        config_data = {
            "sequences": yaml_sequences
        }
        
        # 如果启用亲和力预测，添加properties部分
        if enable_affinity_prediction and affinity_binder_id.strip():
            config_data["properties"] = [
                {
                    "affinity": {
                        "binder": affinity_binder_id.strip()
                    }
                }
            ]
        
        with open(yaml_path, 'w') as f:
            yaml.dump(config_data, f, sort_keys=False)

        yield f"✅ YAML 配置文件已生成于: {yaml_path}\n配置了 {len(yaml_sequences)} 个分子链\n", initial_3d_html, "等待中...", "等待中...", None, None, None

        # 3. 构建并运行 boltz 命令
        cmd = [
            "boltz", "predict", str(input_dir),
            "--out_dir", str(output_dir),
            "--recycling_steps", str(recycling_steps),
            "--diffusion_samples", str(diffusion_samples),
            "--output_format", "mmcif", # 使用mmCIF以获得最佳兼容性
            "--override" # 允许覆盖旧结果（在临时目录中通常不需要）
        ]
        
        # 添加GPU配置
        if gpu_count > 1:
            cmd.extend(["--devices", str(gpu_count)])
        
        if use_msa_server:
            cmd.append("--use_msa_server")
        if use_potentials:
            cmd.append("--use_potentials")
            
        gpu_info = f"使用 {gpu_count} 个GPU" if gpu_count > 1 else "使用单GPU"
        yield f"⚙️ 准备运行 Boltz ({gpu_info})...\n命令: {' '.join(cmd)}\n\n", initial_3d_html, "等待中...", "等待中...", None, None, None

        # 使用 Popen 实时流式传输输出
        process = subprocess.Popen(
            cmd, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT, 
            text=True, 
            encoding='utf-8'
        )
        
        log_output = ""
        while True:
            line = process.stdout.readline()
            if not line:
                break
            log_output += line
            yield log_output, initial_3d_html, "运行中...", "运行中...", None, None, None
        
        process.wait()

        if process.returncode != 0:
            final_log = log_output + f"\n\n❌ Boltz 进程以错误码 {process.returncode} 结束。"
            yield final_log, initial_3d_html, "错误", "错误", None, None, None
            return

        final_log = log_output + "\n\n✅ Boltz 预测完成！"
        yield final_log, initial_3d_html, "处理结果中...", "处理结果中...", None, None, None

        # 4. 处理输出文件
        prediction_folder = output_dir / "boltz_results_input/predictions" / config_name
        
        # 找到排名第一的结构文件
        best_structure_file = prediction_folder / f"{config_name}_model_0.cif"
        
        # 找到对应的置信度和亲和力文件
        confidence_file = prediction_folder / f"confidence_{config_name}_model_0.json"
        affinity_file = prediction_folder / f"affinity_{config_name}.json"

        structure_html_content = "<div style='height: 600px; display: flex; align-items: center; justify-content: center;'><p>未找到结构文件。</p></div>"
        best_structure_file_path_for_download = None

        if best_structure_file.exists():
            try:
                mmcif_bytes = best_structure_file.read_bytes()
                mmcif_base64 = base64.b64encode(mmcif_bytes).decode('utf-8')
                structure_html_content = get_molstar_html(mmcif_base64)
                best_structure_file_path_for_download = str(best_structure_file)
            except Exception as e:
                final_log += f"\n\n❌ 错误：处理结构文件以在Mol*中显示时出错: {e}"
                structure_html_content = f"<div style='height: 600px; display: flex; align-items: center; justify-content: center;'><p>处理结构文件时出错: {e}</p></div>"
        else:
            final_log += "\n\n❌ 错误：找不到预测的结构文件！"
            structure_html_content = get_molstar_html("") # 显示空的Mol*查看器
            
        # 读取JSON数据
        confidence_data = {}
        if confidence_file.exists():
            with open(confidence_file, 'r') as f:
                confidence_data = json.load(f)
        
        affinity_data = {}
        if affinity_file.exists():
            with open(affinity_file, 'r') as f:
                affinity_data = json.load(f)
        
        # 格式化Markdown输出
        confidence_md = create_formatted_confidence_markdown(confidence_data)
        affinity_md = create_formatted_affinity_markdown(affinity_data)
        
        yield (final_log + "\n\n🎉 结果已加载。", 
               structure_html_content, 
               confidence_md, 
               affinity_md,
               best_structure_file_path_for_download,
               str(confidence_file) if confidence_file.exists() else None,
               str(affinity_file) if affinity_file.exists() else None
              )

    except FileNotFoundError:
        yield ("❌ 错误: `boltz` 命令未找到。\n"
               "请确保您已经安装了 `boltz-prediction`并且 `boltz` 在您的系统PATH中。", 
               initial_3d_html, "错误", "错误", None, None, None)
    except Exception as e:
        yield f"❌ 发生意外错误: {e}", initial_3d_html, "错误", "错误", None, None, None
    # 注意：我们不清理临时目录，因为 Gradio 需要从那里提供文件下载。
    # Gradio 会在会话结束后自动处理临时文件。

# --- Gradio 界面 ---
with gr.Blocks(theme=gr.themes.Base()) as demo:
    gr.Markdown(
        """
        # Boltz 生物分子复合物预测工具
        支持蛋白质、DNA、RNA和配体的多条链结构预测，包括结合亲和力计算。
        
        **支持的分子类型：**
        - 🧬 **蛋白质** - 氨基酸序列 (标准单字母代码)
        - 🧬 **DNA** - 脱氧核苷酸序列 (A, T, C, G)
        - 🧬 **RNA** - 核糖核苷酸序列 (A, U, C, G)
        - 💊 **配体** - SMILES字符串或CCD代码 (如ATP, GTP, SAH等)
        
        **功能特色：**
        - ✅ 支持多条链复合物预测
        - ✅ 自动生成MSA (多序列比对)
        - ✅ 结构置信度评估
        - ✅ 结合亲和力预测 (适用于配体)
        - ✅ 交互式3D结构查看器
        """
    )
    
    # 添加帮助信息折叠面板
    with gr.Accordion("📖 使用说明", open=False):
        gr.Markdown(
            """
            ### 基本使用流程：
            1. **配置分子序列**：
               - 手动添加：输入链ID、选择分子类型、输入序列，然后点击"添加分子"
               - 快速开始：点击示例按钮快速加载预设配置
               - 删除配置：在"删除指定链ID"框中输入要删除的链ID（支持用逗号分隔多个，如 A,B,C），然后点击"删除"
               - 清空重置：点击"清空所有"按钮清除所有已配置的分子
               
            2. **选择预测选项**：
               - 建议保持"使用在线MSA服务器"选项开启
               - 可选择启用"推理势能"以提高物理真实性
               - 如需计算配体结合亲和力，勾选相应选项并指定配体链ID
               
            3. **开始预测**：点击"开始预测"按钮，等待计算完成
            
            4. **查看结果**：在不同标签页中查看运行日志、3D结构、置信度和亲和力分数
            
            ### 输入格式说明：
            - **蛋白质序列**：使用标准氨基酸单字母代码，如 `MKITIGSGVSAAKKFV...`
            - **DNA序列**：使用核苷酸字母，如 `ATCGATCGATCG`
            - **RNA序列**：使用核苷酸字母，如 `AUCGAUCGAUCG`
            - **SMILES配体**：化学结构的SMILES表示，如 `C1=CC=C(C=C1)C(=O)O` (苯甲酸)
            - **CCD配体**：PDB化学组分字典代码，如 `ATP`, `GTP`, `SAH`
            
            ### 注意事项：
            - 每个链ID必须唯一
            - 建议使用简短易识别的链ID (如A, B, C, L1, L2等)
            - 亲和力预测仅支持配体分子作为结合物
            - 复杂结构的预测时间较长，请耐心等待
            - GPU数量设置：更多GPU可以加速预测，但也会消耗更多显存
            - 临时文件自动保存在当前目录的tmp文件夹中，便于查看和调试
            """
        )

    with gr.Row():
        with gr.Column(scale=1):
            gr.Markdown("### 步骤 1: 配置分子序列")
            
            # 添加示例按钮
            with gr.Row():
                example1_btn = gr.Button("📝 示例1: 蛋白质-配体", variant="secondary", size="sm")
                example2_btn = gr.Button("📝 示例2: 蛋白质-DNA", variant="secondary", size="sm")
                example3_btn = gr.Button("📝 示例3: 多链复合物", variant="secondary", size="sm")
            
            # 用于存储序列配置的状态
            sequences_state = gr.State([])
            
            with gr.Row():
                chain_id_input = gr.Textbox(label="链ID", placeholder="例如: A, B, C, L1", scale=1)
                mol_type_input = gr.Radio(
                    ["蛋白质", "DNA", "RNA", "配体(SMILES)", "配体(CCD)"], 
                    label="分子类型", 
                    value="蛋白质",
                    scale=1
                )
            
            sequence_input = gr.Textbox(
                label="序列/标识符", 
                placeholder="输入氨基酸序列、核苷酸序列、SMILES字符串或CCD代码...",
                lines=3
            )
            
            with gr.Row():
                add_sequence_btn = gr.Button("➕ 添加分子", variant="secondary")
                clear_sequences_btn = gr.Button("🗑️ 清空所有", variant="secondary")
            
            with gr.Row():
                delete_chain_id = gr.Textbox(
                    label="删除指定链ID", 
                    placeholder="输入要删除的链ID，多个用逗号分隔 (例如: A,B,C)",
                    scale=3
                )
                delete_specific_btn = gr.Button("❌ 删除", variant="secondary", scale=1)
            
            # 显示当前配置的序列
            sequences_display = gr.Dataframe(
                headers=["链ID", "分子类型", "序列/标识符"],
                datatype=["str", "str", "str"],
                label="已配置的分子序列",
                interactive=False
            )

            gr.Markdown("### 步骤 2: 配置预测选项")
            use_msa_server = gr.Checkbox(label="使用在线 MSA 服务器 (推荐)", value=True)
            use_potentials = gr.Checkbox(label="使用推理势能 (提高物理真实性)", value=False)
            
            # 亲和力预测选项
            with gr.Group():
                gr.Markdown("#### 亲和力预测 (可选)")
                enable_affinity = gr.Checkbox(label="启用结合亲和力预测", value=False)
                affinity_binder_id = gr.Textbox(
                    label="结合分子链ID", 
                    placeholder="例如: L1 (必须是上面已添加的链ID)",
                    visible=False
                )
            
            # 高级选项
            with gr.Accordion("高级选项", open=False):
                # GPU配置
                available_gpus = get_available_gpus()
                gpu_count = gr.Slider(
                    minimum=1, maximum=available_gpus, value=available_gpus, step=1,
                    label=f"使用GPU数量 (检测到 {available_gpus} 个GPU)",
                    info=f"选择用于预测的GPU数量。默认使用所有可用GPU ({available_gpus}个)"
                )
                
                recycling_steps = gr.Slider(
                    minimum=1, maximum=10, value=3, step=1, 
                    label="循环步数 (Recycling Steps)",
                    info="更多的步数可能提高精度，但会增加时间。默认: 3"
                )
                diffusion_samples = gr.Slider(
                    minimum=1, maximum=10, value=1, step=1,
                    label="扩散样本数 (Diffusion Samples)",
                    info="生成多个候选结构以选择最优。默认: 1"
                )
            
            run_button = gr.Button("🚀 开始预测", variant="primary", size="lg")

        with gr.Column(scale=2):
            gr.Markdown("### 步骤 3: 查看结果")
            
            with gr.Tabs():
                with gr.TabItem("📈 运行日志"):
                    status_log = gr.Textbox(label="状态和日志", lines=15, interactive=False)
                with gr.TabItem("🔬 3D 结构"):
                    model_3d_view = gr.HTML(
                        label="最佳预测结构 (排名 1)", 
                        value="<div style='height: 600px; display: flex; align-items: center; justify-content: center;'><p>等待预测开始...</p></div>"
                    )
                with gr.TabItem("📊 置信度分数"):
                    confidence_output = gr.Markdown("预测完成后，此处将显示置信度分数。")
                with gr.TabItem("💞 亲和力分数"):
                    affinity_output = gr.Markdown("预测完成后，此处将显示亲和力分数。")

            gr.Markdown("#### 📂 下载结果文件")
            with gr.Row():
                download_structure = gr.File(label="下载结构 (.cif)")
                download_confidence = gr.File(label="下载置信度 (.json)")
                download_affinity = gr.File(label="下载亲和力 (.json)")

    # 定义用于管理序列状态的辅助函数
    def add_sequence(chain_id, mol_type, sequence, current_sequences):
        """添加新的序列配置到状态中"""
        if not chain_id.strip() or not sequence.strip():
            return current_sequences, current_sequences, "错误：链ID和序列不能为空"
        
        # 检查链ID是否已存在
        for seq in current_sequences:
            if seq["chain_id"] == chain_id.strip():
                return current_sequences, current_sequences, f"错误：链ID '{chain_id}' 已存在"
        
        new_sequence = {
            "chain_id": chain_id.strip(),
            "mol_type": mol_type,
            "sequence": sequence.strip()
        }
        updated_sequences = current_sequences + [new_sequence]
        
        # 转换为显示格式
        display_data = [[seq["chain_id"], seq["mol_type"], seq["sequence"][:50] + "..." if len(seq["sequence"]) > 50 else seq["sequence"]] 
                       for seq in updated_sequences]
        
        return updated_sequences, display_data, f"✅ 成功添加分子：{chain_id} ({mol_type})"
    
    def clear_sequences():
        """清空所有序列配置"""
        return [], [], "✅ 已清空所有分子配置"
    
    def delete_specific_sequence(current_sequences, chain_ids_to_delete):
        """删除指定链ID的序列配置，支持多个链ID用逗号分隔"""
        if not chain_ids_to_delete.strip():
            return current_sequences, current_sequences, "请输入要删除的链ID"
        
        # 解析链ID列表
        chain_ids = [cid.strip() for cid in chain_ids_to_delete.split(',') if cid.strip()]
        
        if not chain_ids:
            return current_sequences, current_sequences, "请输入有效的链ID"
        
        # 查找要删除的序列
        remaining_sequences = [seq for seq in current_sequences if seq["chain_id"] not in chain_ids]
        
        deleted_count = len(current_sequences) - len(remaining_sequences)
        
        if deleted_count == 0:
            return current_sequences, current_sequences, f"错误：未找到指定的链ID: {', '.join(chain_ids)}"
        
        # 转换为显示格式
        display_data = [[seq["chain_id"], seq["mol_type"], seq["sequence"][:50] + "..." if len(seq["sequence"]) > 50 else seq["sequence"]] 
                       for seq in remaining_sequences]
        
        if deleted_count == 1:
            return remaining_sequences, display_data, f"✅ 成功删除链ID '{chain_ids[0]}' 的分子配置"
        else:
            return remaining_sequences, display_data, f"✅ 成功删除 {deleted_count} 个分子配置: {', '.join(chain_ids[:deleted_count])}"
    
    def toggle_affinity_options(enable_affinity):
        """切换亲和力预测选项的可见性"""
        return gr.update(visible=enable_affinity)
    
    def load_example1():
        """加载示例1：蛋白质-配体复合物"""
        example_sequences = [
            {"chain_id": "A", "mol_type": "蛋白质", "sequence": "MKITIGSGVSAAKKFVGLKQPGRYDYKVLAYPIAVEALSLIYNKDLLPNPPKTWEEIPALDKELKAFDISTEELSA"},
            {"chain_id": "L", "mol_type": "配体(SMILES)", "sequence": "C1=CC=C(C=C1)C(=O)O"}
        ]
        display_data = [[seq["chain_id"], seq["mol_type"], seq["sequence"][:50] + "..." if len(seq["sequence"]) > 50 else seq["sequence"]] 
                       for seq in example_sequences]
        return example_sequences, display_data, "✅ 已加载示例1：蛋白质-配体复合物"
    
    def load_example2():
        """加载示例2：蛋白质-DNA复合物"""
        example_sequences = [
            {"chain_id": "A", "mol_type": "蛋白质", "sequence": "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPGCMSKCVLS"},
            {"chain_id": "D", "mol_type": "DNA", "sequence": "ATCGATCGATCGATCG"}
        ]
        display_data = [[seq["chain_id"], seq["mol_type"], seq["sequence"][:50] + "..." if len(seq["sequence"]) > 50 else seq["sequence"]] 
                       for seq in example_sequences]
        return example_sequences, display_data, "✅ 已加载示例2：蛋白质-DNA复合物"
    
    def load_example3():
        """加载示例3：多链复合物 (蛋白质二聚体 + RNA + 配体)"""
        example_sequences = [
            {"chain_id": "A", "mol_type": "蛋白质", "sequence": "MKITIGSGVSAAKKFVGLKQPGRYDYKVLAYPIAVEALSLIYNKDLLPNPPKTWEEIPALDKELKAFDISTEELSA"},
            {"chain_id": "B", "mol_type": "蛋白质", "sequence": "MKITIGSGVSAAKKFVGLKQPGRYDYKVLAYPIAVEALSLIYNKDLLPNPPKTWEEIPALDKELKAFDISTEELSA"},
            {"chain_id": "R", "mol_type": "RNA", "sequence": "AUCGAUCGAUCGAUCG"},
            {"chain_id": "L1", "mol_type": "配体(CCD)", "sequence": "ATP"}
        ]
        display_data = [[seq["chain_id"], seq["mol_type"], seq["sequence"][:50] + "..." if len(seq["sequence"]) > 50 else seq["sequence"]] 
                       for seq in example_sequences]
        return example_sequences, display_data, "✅ 已加载示例3：多链复合物 (蛋白质二聚体 + RNA + ATP)"

    # 事件绑定
    # 示例按钮事件
    example1_btn.click(
        fn=load_example1,
        outputs=[sequences_state, sequences_display, status_log]
    )
    
    example2_btn.click(
        fn=load_example2,
        outputs=[sequences_state, sequences_display, status_log]
    )
    
    example3_btn.click(
        fn=load_example3,
        outputs=[sequences_state, sequences_display, status_log]
    )
    
    # 序列管理事件
    add_sequence_btn.click(
        fn=add_sequence,
        inputs=[chain_id_input, mol_type_input, sequence_input, sequences_state],
        outputs=[sequences_state, sequences_display, status_log]
    )
    
    clear_sequences_btn.click(
        fn=clear_sequences,
        outputs=[sequences_state, sequences_display, status_log]
    )
    
    delete_specific_btn.click(
        fn=delete_specific_sequence,
        inputs=[sequences_state, delete_chain_id],
        outputs=[sequences_state, sequences_display, status_log]
    )
    
    enable_affinity.change(
        fn=toggle_affinity_options,
        inputs=[enable_affinity],
        outputs=[affinity_binder_id]
    )

    # 将预测按钮点击事件连接到处理函数
    run_button.click(
        fn=run_boltz_prediction,
        inputs=[
            sequences_state,
            use_msa_server,
            use_potentials,
            recycling_steps, 
            diffusion_samples,
            enable_affinity,
            affinity_binder_id,
            gpu_count
        ],
        outputs=[
            status_log,
            model_3d_view,
            confidence_output,
            affinity_output,
            download_structure,
            download_confidence,
            download_affinity
        ]
    )

if __name__ == "__main__":
    demo.launch(debug=True, server_name="0.0.0.0")