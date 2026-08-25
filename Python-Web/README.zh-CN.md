# MagneticKP Python/Web 版

简体中文 | [English](README.md)

本目录包含 MagneticKP 的 Python/Web 稳定发行版：

* `magnetickp/`：Python 软件包和中英双语 Web 界面；
* `tests/`：Python 与浏览器端测试；
* `setup.py` 和 `MANIFEST.in`：打包配置。

MagneticKP 和随附的 Python MSGCorep 发行包都是纯 Python 软件包，并包含各自的浏览器资源和数据库。它们的依赖项包括与平台相关的 NumPy、SciPy 和 spglib wheel，但安装后的应用不需要 Mathematica。

`dist/` 目录包含相互匹配的 MagneticKP 和 MSGCorep wheel。建议先创建并激活虚拟环境，避免修改系统 Python 环境：

```console
cd /你的路径/MagneticKP1.2/MagneticKP/Python-Web
python3 -m venv env
source env/bin/activate
python -m pip install --upgrade pip
python -m pip install --find-links dist \
  dist/magnetickp-1.2.0-py3-none-any.whl
magnetickp-web
```

上述安装命令会一次安装 Python 核心、MSGCorep 目录、Web 应用、静态绘图和 VASP 拟合。使用完毕后运行 `deactivate` 退出虚拟环境。以后再次使用 MagneticKP 时，先进入本目录并运行 `source env/bin/activate` 即可。

应用默认监听 `http://127.0.0.1:8000`。如需修改监听地址或端口，请在启动前设置 `MAGNETICKP_WEB_HOST` 或 `MAGNETICKP_WEB_PORT`。

构造好 `ham` 和 `path` 后，只需一次调用：

```python
from magnetickp import bandplot
bandplot(ham, path, parameter_values)
```

构建可重新分发的源码包和 wheel：

```console
python3 -m pip install build
python3 -m build
```

分发完整的本地安装包之前，请将 `msgcorep-0.1.0a0-py3-none-any.whl` 与生成的 MagneticKP wheel 放在同一个 `dist/` 目录中。

开发时，从同一本地 wheel 目录安装可编辑软件包和测试依赖，然后运行测试：

```console
python3 -m pip install --find-links dist -e '.[test]'
python3 -m pytest -q
```
