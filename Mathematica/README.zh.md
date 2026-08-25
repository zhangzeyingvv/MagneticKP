# MagneticKP for Wolfram Language

[English](README.md) | [简体中文](README.zh.md)

本目录包含 MagneticKP 的 Mathematica 稳定发行版。可直接安装的程序包位于
[`build/MagneticKP-1.2.0.paclet`](build/MagneticKP-1.2.0.paclet)。

## 安装 Paclet

打开 Mathematica，运行以下代码，然后选择已下载的
`MagneticKP-1.2.0.paclet` 文件：

```wl
pacletFile = SystemDialogInput["FileOpen"];
PacletInstall[pacletFile, ForceVersionInstall -> True]
```

`ForceVersionInstall -> True` 会替换已经安装的同版本程序包，适合重新安装发行包。
安装完成后请启动一个新的内核，再载入 MagneticKP：

```wl
Needs["MagneticKP`"]
```

可以用以下代码确认安装结果并打开 MagneticKP 的主帮助页：

```wl
PacletFind["MagneticKP"]
Documentation`HelpLookup["paclet:MagneticKP/guide/MagneticKP"]
```

安装后的帮助系统包含英文和简体中文页面。函数页面和教程中保存了示例的运行
结果，用户无需先运行代码即可查看预期输出。不需要手动将文件复制到
`$UserBaseDirectory/Applications`。

## 开发与测试

其余目录包含源代码、示例、文档和验证工具：

* `MagneticKP/`：Paclet 源文件；
* `Tests/`：单元测试、回归测试、集成测试和穷举测试；
* `Benchmarks/`：基准测试与原型脚本；
* `Examples/`：示例笔记本和样例数据；
* `Documentation/`：PDF 手册源文件及生成的手册；
* `Developer/Documentation/`：帮助文档生成规范、翻译、预期文本以及笔记本构建脚本。

## 生成帮助笔记本

`MagneticKP/Documentation/` 下的文件是自动生成的帮助笔记本，不要直接编辑。
请在维护工作区的根目录运行以下命令，依次生成英文规范、简体中文翻译和规范化
笔记本：

```console
env OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /Applications/Mathematica.app/Contents/MacOS/WolframKernel \
  -noinit -script MagneticKP/Mathematica/Developer/Documentation/GenerateDocumentationSources.wls
env OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /Applications/Mathematica.app/Contents/MacOS/WolframKernel \
  -noinit -script MagneticKP/Mathematica/Developer/Documentation/GenerateChineseDocumentationSources.wls
MagneticKP/Mathematica/Developer/Documentation/BuildNotebooks.sh
```

英文内容和页面结构在 `GenerateDocumentationSources.wls` 中编辑，简体中文文字在
`ChineseTranslations.wl` 中编辑。`BuildNotebooks.sh` 会规范化并验证所有生成的
笔记本。

## 构建 Paclet

在维护工作区的根目录运行以下命令构建可分发的 Paclet：

```console
wolframscript -file MagneticKP/Mathematica/BuildPaclet.wls
```

生成的安装包位于 `MagneticKP/Mathematica/build/`。

构建完成后，运行以下命令检查双语页面覆盖、帮助链接、生成笔记本的解析、Paclet
注册以及文档中的核心示例：

```console
env OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /Applications/Mathematica.app/Contents/MacOS/WolframKernel \
  -noinit -script MagneticKP/Mathematica/Developer/Documentation/VerifyDocumentation.wls
```

在维护工作区的根目录运行常规测试：

```console
wolframscript -file MagneticKP/Mathematica/Tests/RunTests.wls
```
