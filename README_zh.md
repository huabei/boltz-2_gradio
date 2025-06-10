## Introduction to Boltz

Boltz 是一个用于生物分子相互作用预测的模型家族。Boltz-1 是第一个接近 AlphaFold3 精度的完全开源模型。最新的 Boltz-2 是一个新的生物分子基础模型，它通过联合建模复杂的结构和结合亲和力，超越了 AlphaFold3 和 Boltz-1，这是实现精确分子设计的关键组成部分。Boltz-2 是第一个在精度上接近基于物理的自由能微扰 (FEP) 方法的深度学习模型，同时运行速度提高了 1000 倍——使得精确的计算机筛选在早期药物发现中变得实用。

所有的代码和权重都在 MIT 许可下提供，可免费用于学术和商业用途。更多关于模型的信息，请参阅 [Boltz-1](https://doi.org/10.1101/2024.11.19.624167) 和 [Boltz-2](https://bit.ly/boltz2-pdf) 技术报告。

## Boltz-2 Gradio Interface

本项目提供了一个基于 [Gradio](https://www.gradio.app/) 的用户友好界面，用于与 Boltz-2 模型进行交互。通过这个界面，用户可以更方便地输入生物分子信息并获取预测结果，如结合亲和力等，而无需直接操作命令行。

# Usage

## 安装

我们建议在一个新的 Python 环境中安装。

1.  **安装 Boltz-2：**
    请遵循 [Boltz 官方仓库](https://github.com/jwohlwend/boltz) 的指引安装 Boltz-2。通常可以通过 pip 安装：
    ```bash
    pip install boltz -U
    ```
    或者从 GitHub 安装以获取最新更新：
    ```bash
    git clone https://github.com/jwohlwend/boltz.git
    cd boltz; pip install -e .
    ```

2.  **安装 Gradio：**
    ```bash
    pip install gradio
    ```

3.  **克隆本项目：**
    ```bash
    git clone https://github.com/huabei/boltz-2_gradio.git
    cd boltz-2_gradio
    ```

## 运行应用

完成安装后，在项目根目录下运行 `app.py` 文件：

```bash
python app.py
```

之后，在浏览器中打开输出的 Gradio 链接 (通常是 `http://127.0.0.1:7860`) 即可开始使用。