import gradio as gr
import yaml
import subprocess
import tempfile
import shutil
import os
import json
from pathlib import Path
import base64 # 新增导入

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
    protein_sequence, 
    ligand_type, 
    ligand_identifier, 
    use_msa_server,
    use_potentials,
    recycling_steps, 
    diffusion_samples
):
    """
    一个完整的函数，用于生成YAML，运行Boltz，并处理输出。
    """
    # 1. 输入验证
    if not protein_sequence.strip():
        return "错误：蛋白质序列不能为空。", None, None, None, None, None, None
    if not ligand_identifier.strip():
        return f"错误：{ligand_type} 不能为空。", None, None, None, None, None, None

    # 创建一个临时目录来存放所有文件
    run_dir = tempfile.mkdtemp(prefix="boltz_gradio_run_")
    # run_dir = 'tmp_boltz_run'  # 使用固定目录以便于调试和查看结果
    # Path(run_dir).mkdir(exist_ok=True)
    
    initial_3d_html = "<div style='height: 600px; display: flex; align-items: center; justify-content: center;'><p>等待预测结束...</p></div>"

    try:
        # 2. 生成 YAML 配置文件
        input_dir = Path(run_dir) / "input"
        output_dir = Path(run_dir) / "output"
        input_dir.mkdir()
        output_dir.mkdir()

        config_name = "prediction_config"
        yaml_path = input_dir / f"{config_name}.yaml"
        
        # 定义配体部分
        ligand_entry = {"id": "L"} # 使用固定的链ID 'L'
        if ligand_type == "SMILES":
            ligand_entry["smiles"] = ligand_identifier
        else: # CCD Code
            ligand_entry["ccd"] = ligand_identifier

        # 构建完整的YAML结构
        config_data = {
            "sequences": [
                {
                    "protein": {
                        "id": "A", # 使用固定的链ID 'A'
                        "sequence": protein_sequence.strip()
                    }
                },
                {
                    "ligand": ligand_entry
                }
            ],
            "properties": [
                {
                    "affinity": {
                        "binder": "L" # 指定配体链ID以计算亲和力
                    }
                }
            ]
        }
        
        with open(yaml_path, 'w') as f:
            yaml.dump(config_data, f, sort_keys=False)

        yield f"✅ YAML 配置文件已生成于: {yaml_path}\n", initial_3d_html, "等待中...", "等待中...", None, None, None

        # 3. 构建并运行 boltz 命令
        cmd = [
            "boltz", "predict", str(input_dir),
            "--out_dir", str(output_dir),
            "--recycling_steps", str(recycling_steps),
            "--diffusion_samples", str(diffusion_samples),
            "--output_format", "mmcif", # 使用mmCIF以获得最佳兼容性
            "--override" # 允许覆盖旧结果（在临时目录中通常不需要）
        ]
        if use_msa_server:
            cmd.append("--use_msa_server")
        if use_potentials:
            cmd.append("--use_potentials")
            
        yield f"⚙️ 准备运行 Boltz...\n命令: {' '.join(cmd)}\n\n", initial_3d_html, "等待中...", "等待中...", None, None, None

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
        # Boltz 蛋白质-配体复合物预测工具
        一个简单的界面，用于预测蛋白质-配体复合物的结构和结合亲和力。
        """
    )

    with gr.Row():
        with gr.Column(scale=1):
            gr.Markdown("### 步骤 1: 输入分子信息")
            protein_seq = gr.Textbox(
                label="蛋白质序列 (链 A)", 
                placeholder="输入单条蛋白质链的氨基酸序列...",
                lines=5
            )
            ligand_type = gr.Radio(
                ["SMILES", "CCD Code"], 
                label="配体标识符类型", 
                value="SMILES"
            )
            ligand_id = gr.Textbox(
                label="配体 SMILES 字符串 (链 L)", 
                placeholder="例如: C1=CC=C(C=C1)C(=O)O"
            )

            def update_ligand_label(choice):
                if choice == "SMILES":
                    return gr.Textbox(label="配体 SMILES 字符串 (链 L)", placeholder="例如: C1=CC=C(C=C1)C(=O)O")
                else:
                    return gr.Textbox(label="配体 CCD Code (链 L)", placeholder="例如: ATP, SAH")

            ligand_type.change(fn=update_ligand_label, inputs=ligand_type, outputs=ligand_id)

            gr.Markdown("### 步骤 2: 配置预测选项")
            use_msa_server = gr.Checkbox(label="使用在线 MSA 服务器 (必须)", value=True)
            use_potentials = gr.Checkbox(label="使用推理势能 (提高物理真实性)", value=False)
            
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
            
            run_button = gr.Button("🚀 开始预测", variant="primary")

        with gr.Column(scale=2):
            gr.Markdown("### 步骤 3: 查看结果")
            
            with gr.Tabs():
                with gr.TabItem("📈 运行日志"):
                    status_log = gr.Textbox(label="状态和日志", lines=15, interactive=False)
                with gr.TabItem("🔬 3D 结构"):
                    model_3d_view = gr.HTML(
                        label="最佳预测结构 (排名 1)", 
                        # value=get_initial_molstar_html() # 初始加载默认 example.cif 或空查看器
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

    # 将按钮点击事件连接到处理函数
    run_button.click(
        fn=run_boltz_prediction,
        inputs=[
            protein_seq, 
            ligand_type, 
            ligand_id, 
            use_msa_server,
            use_potentials,
            recycling_steps, 
            diffusion_samples
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