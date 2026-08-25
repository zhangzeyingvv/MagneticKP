# MagneticKP

[English](README.md)

MagneticKP 用于根据对称操作和（共）表示矩阵构造 k·p 哈密顿量，提供 Mathematica、Python 和 Web 三种使用方式。

## 安装方法

### Mathematica

从项目的 [Releases 页面](https://github.com/zhangzeyingvv/MagneticKP/releases) 下载 `MagneticKP-1.2.0.paclet`，然后在 Mathematica 中运行：

```wl
PacletInstall[SystemDialogInput["FileOpen"]]
Needs["MagneticKP`"]
```

在弹出的窗口中选择下载的 `.paclet` 文件即可。需要 Mathematica 12.1 或更高版本。

### Python/Web

进入 `Python-Web` 目录。建议在该目录中创建并使用虚拟环境：

```bash
python3 -m venv env
source env/bin/activate
python -m pip install --upgrade pip
python -m pip install --find-links dist dist/magnetickp-1.2.0-py3-none-any.whl
magnetickp-web
```

安装命令会从同一个 `dist/` 目录自动安装 `msgcorep`，无需单独安装。然后在浏览器打开 `http://127.0.0.1:8000`。

使用完毕后，运行以下命令退出虚拟环境：

```bash
deactivate
```

以后再次使用 MagneticKP 时，进入 `Python-Web` 目录，先运行 `source env/bin/activate`，再启动 Python 或 `magnetickp-web`。

## MagneticKP 的功能

* 构造任意阶 k·p 哈密顿量；
* 支持磁性、非磁性和自旋轨道耦合体系；
* 支持内部连续自旋旋转 `CInfinity`；
* 通过 MSGCorep 选择磁空间群、k 点和共表示；
* 根据所选 k 点的小陪群自动计算局域高对称 k 路径；
* 绘制能带，并结合第一性原理能带拟合参数；
* 提供 Mathematica、Python 和 Web 界面。

## 如何打开说明书

安装 Paclet 后，在 Mathematica 中运行以下代码打开内置帮助：

```wl
Documentation`HelpLookup["paclet:MagneticKP/guide/MagneticKP"]
```

内置帮助包含英文和简体中文页面。

Python/Web 版请运行 `magnetickp-web`，然后用浏览器打开 `http://127.0.0.1:8000/help/`。

更多详细说明请查看 [Mathematica/README.md](Mathematica/README.md) 和 [Python-Web/README.md](Python-Web/README.md)。

## 简单例子

### Mathematica

```wl
input = <|
  "Unitary" -> <|
    "C2z" -> {PauliMatrix[3], {-kx, -ky, kz}}
  |>
|>;

kpHam[2, input]
```

### Python

```python
from sympy import Matrix
from magnetickp import kpHam

input_data = {
    "Unitary": {
        "C2z": (
            Matrix([[1, 0], [0, -1]]),
            Matrix([[-1, 0, 0], [0, -1, 0], [0, 0, 1]]),
        )
    }
}

result = kpHam(input_data).getkpHam(2)
print(result["ham"])
```

## 更新记录

### v1.2.0 — 2026/8/20

* 新增 Python/Web 磁空间群和共表示选择界面；
* 新增根据所选 k 点小陪群自动计算局域高对称 k 路径；
* Mathematica和python+web增加了大量帮助文档；

### v1.1.0 — 2025/10/31

* 提高动量变换矩阵的计算速度；
* 新增 Reynolds 求解方式；
* 支持连续自旋旋转 `CInfinity`；

### v1.04 — 2023/10/31

* 新增 `vaspEig`、`bandManipulateWithVASP` 和 `fittingKP`，用于拟合 k·p 哈密顿量。

### v1.03 — 2023/03/26

* 新增 Python 版本。

### v1.02 — 2022/10/15

* 新增 `bandManipulate` 和 `bandplot` 能带绘图功能；
* 支持 MSGCorep v1.0.0。

### v1.01 — 2022/06/30

* 支持一维和二维 k·p 哈密顿量。
