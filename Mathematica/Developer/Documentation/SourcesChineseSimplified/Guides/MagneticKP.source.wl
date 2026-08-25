Notebook[{Cell["MagneticKP", "GuideTitle"],
  Cell["MagneticKP 用于构造满足对称性要求的 k.p 哈密顿量，可连接表示数据库软件包、绘制模型能带，并根据 VASP \
数据拟合模型参数。", "GuideAbstract"],
  Cell[CellGroupData[{Cell["", "GuideFunctionsSection"],
     Cell[CellGroupData[{Cell["构造 k.p 模型", "GuideFunctionsSubsection"],
        Cell[TextData[{ButtonBox["kpHam", BaseStyle -> "Link", ButtonData -> "paclet:MagneticKP/ref/kpHam"],
           " — construct a symmetry-allowed k.p Hamiltonian."}], "GuideText"],
        Cell[TextData[{ButtonBox["interfaceRep", BaseStyle -> "Link", ButtonData ->
             "paclet:MagneticKP/ref/interfaceRep"], " — prepare kpHam input from SpaceGroupIrep or MSGCorep."}],
         "GuideText"], Cell[TextData[{ButtonBox["IterativelySimplify", BaseStyle -> "Link",
            ButtonData -> "paclet:MagneticKP/ref/IterativelySimplify"],
           " — call the default construction method directly."}], "GuideText"],
        Cell[TextData[{ButtonBox["GETkpHam", BaseStyle -> "Link", ButtonData -> "paclet:MagneticKP/ref/GETkpHam"],
           " — use the legacy homogeneous-order interface."}], "GuideText"]}, Open]], Cell["\t", "GuideDelimiter"],
     Cell[CellGroupData[{Cell["能带与参数拟合", "GuideFunctionsSubsection"],
        Cell[TextData[{ButtonBox["highSymmetryKPath", BaseStyle -> "Link", ButtonData ->
             "paclet:MagneticKP/ref/highSymmetryKPath"],
           " — generate a local momentum path from little-co-group subgroups."}], "GuideText"],
        Cell[TextData[{ButtonBox["bandManipulate", BaseStyle -> "Link", ButtonData ->
             "paclet:MagneticKP/ref/bandManipulate"], " — explore a model with interactive parameter controls."}],
         "GuideText"], Cell[TextData[{ButtonBox["bandplot", BaseStyle -> "Link", ButtonData ->
             "paclet:MagneticKP/ref/bandplot"], " — plot a numerical model along a labeled momentum path."}],
         "GuideText"], Cell[TextData[{ButtonBox["vaspEig", BaseStyle -> "Link", ButtonData ->
             "paclet:MagneticKP/ref/vaspEig"], " — import selected bands from a VASP EIGENVAL file."}], "GuideText"],
        Cell[TextData[{ButtonBox["bandManipulateWithVASP", BaseStyle -> "Link", ButtonData ->
             "paclet:MagneticKP/ref/bandManipulateWithVASP"], " — compare model bands with imported VASP data."}],
         "GuideText"], Cell[TextData[{ButtonBox["fittingKP", BaseStyle -> "Link", ButtonData ->
             "paclet:MagneticKP/ref/fittingKP"], " — fit model parameters to selected band data."}], "GuideText"]},
       Open]], Cell["\t", "GuideDelimiter"], Cell[CellGroupData[{Cell["实用函数", "GuideFunctionsSubsection"],
        Cell[TextData[{ButtonBox["getGenerator", BaseStyle -> "Link", ButtonData ->
             "paclet:MagneticKP/ref/getGenerator"], " — find a generating subset of a finite group."}], "GuideText"],
        Cell[TextData[{ButtonBox["directSum", BaseStyle -> "Link", ButtonData -> "paclet:MagneticKP/ref/directSum"],
           " — form a block-diagonal direct sum."}], "GuideText"]}, Open]]}, Open]],
  Cell[CellGroupData[{Cell["教程", "GuideTutorialsSection"],
     Cell[TextData[{ButtonBox["MagneticKP 入门", BaseStyle -> "Link", ButtonData ->
          "paclet:MagneticKP/tutorial/GettingStarted"]}], "GuideTutorial"]}, Open]],
  Cell[CellGroupData[{Cell["历史", "HistorySection"],
     Cell[TextData[{"引入版本：", Cell["1.1.0", "HistoryData", CellTags -> "New"], " | 修改版本：",
        Cell["1.2.0", "HistoryData", CellTags -> "Modified"], " | 废弃版本：",
        Cell[" ", "HistoryData", CellTags -> "Obsolete"]}], "History"]}, Closed]],
  Cell[CellGroupData[{Cell["分类", "CategorizationSection"], Cell["Guide", "Categorization",
      CellLabel -> "Entity Type"], Cell["MagneticKP", "Categorization", CellLabel -> "Paclet Name"],
     Cell["MagneticKP`", "Categorization", CellLabel -> "Context"], Cell["MagneticKP/guide/MagneticKP",
      "Categorization", CellLabel -> "URI"]}, Closed]],
  Cell[CellGroupData[{Cell["关键词", "KeywordsSection"], Cell["k.p 哈密顿量", "Keywords"],
     Cell["对称性", "Keywords"], Cell["band structure", "Keywords"], Cell["VASP", "Keywords"]}, Closed]]},
 TaggingRules -> <|"Paclet" -> "MagneticKP"|>, WindowTitle -> "MagneticKP",
 StyleDefinitions -> FrontEnd`FileName[{"Wolfram"}, "GuidePageStylesExt.nb", CharacterEncoding -> "UTF-8"]]