[中文版](README_zh.md)

## Introduction to Boltz

Boltz is a family of models for biomolecular interaction prediction. Boltz-1 was the first fully open-source model to approach AlphaFold3 accuracy. The latest work, Boltz-2, is a new biomolecular foundation model that goes beyond AlphaFold3 and Boltz-1 by jointly modeling complex structures and binding affinities, a critical component towards accurate molecular design. Boltz-2 is the first deep learning model to approach the accuracy of physics-based free-energy perturbation (FEP) methods, while running 1000x faster—making accurate in silico screening practical for early-stage drug discovery.

All the code and weights are provided under the MIT license, making them freely available for both academic and commercial uses. For more information about the model, see the [Boltz-1](https://doi.org/10.1101/2024.11.19.624167) and [Boltz-2](https://bit.ly/boltz2-pdf) technical reports.

## Boltz-2 Gradio Interface

This project provides a user-friendly interface based on [Gradio](https://www.gradio.app/) for interacting with the Boltz-2 model. Through this interface, users can more conveniently input biomolecular information and obtain prediction results, such as binding affinities, without needing to operate the command line directly.

# Usage

## Installation

We recommend installing in a new Python environment.

1.  **Install Boltz-2:**
    Please follow the instructions in the [Boltz official repository](https://github.com/jwohlwend/boltz) to install Boltz-2. Typically, it can be installed via pip:
    ```bash
    pip install boltz -U
    ```
    Or install from GitHub for the latest updates:
    ```bash
    git clone https://github.com/jwohlwend/boltz.git
    cd boltz; pip install -e .
    ```

2.  **Install Gradio:**
    ```bash
    pip install gradio
    ```

3.  **Clone this project:**
    ```bash
    git clone https://github.com/huabei/boltz-2_gradio.git
    cd boltz-2_gradio
    ```

## Running the Application

After completing the installation, run the `app.py` file in the project root directory:

```bash
python app.py
```

Then, open the Gradio link output in your browser (usually `http://127.0.0.1:7860`) to start using it.