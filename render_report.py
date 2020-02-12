#!/usr/bin/env python3
import sys
import os

sp, html_path, report_type = sys.argv


qc_content = """
            <p>质量控制是生物分析的基本概念之一，用在保证组学测定的数据的重复性和精确性。由于色谱系统与质谱直接与样品接触， 随着分析样品的增多，色谱柱和质谱会逐步的污染，导致信号的漂移，造成测量的系统误差。通过重复使用同一个质控样本(QC样本)来跟踪整个数据采集过程的行为， 已经被大多数的分析化学领域专家推荐和使用。质控样本被用于评估整个质谱数据在采集过程中的信号漂移， 这些漂移进一步能够被精确的算法所识别，校正，提高数据的质量。本流程采用R语言statTarget包的QC-RFSC算法对各个样本的特征（每个代谢物）信号峰进行校正，并记录了每个代谢物的校正效果。QC样本是所有样本取等量混合后的样本，在信号数据采集过程中，在开头，结尾，以及中间部分位点插入QC样本，记录信号漂移情况。所有QC样本都是一样的，如果没有信号漂移，在数据采集过程中，QC样本的信号强度应该是保持不变的。如图3-3所示，校正信号漂移后，如果在PCA图中，QC样本点聚到一起，证明校正效果良好。</p>
        
            <div class="image-block">
                <img src="FiguresTablesForReport/图3-3.png">
                <h6>图3-3 质控样本PCA图</h6>
                <p class="image-description">红色点为校正后的质控样本点（QC样本），蓝色点为检测样本。</p>
            </div>

            <div class="alert alert-primary alert-dismissible">
                <button type="button" class="close" data-dismiss="alert">&times;</button>
                <strong>注：</strong> 数据预处理数据文件：metabolome_pos.xlsx （正离子模式）,metabolome_neg.xlsx （负离子模式）。metabolome_*.xls里面包含“raw”和“normalized”两个工作簿，其中“raw”工作簿是数据预处理前（即解卷积获得的）数据，“normalized”工作簿是经过归一化、QA等预处理后的数据。
            </div>

"""

pathway_nav = """
            <div class="btn-group">
                <button type="button" class="btn btn-info"><a class="buta" href="#a10">十 通路分析</a></button>
                <button type="button" class="btn btn-primary dropdown-toggle dropdown-toggle-split" data-toggle="dropdown">
                    <span class="caret"></span>
                </button>
                <div class="dropdown-menu">
                    <a class="dropdown-item" href="#a10.1">10.1 富集分析</a>
                    <a class="dropdown-item" href="#a10.2">10.2 拓扑分析</a>
                    <a class="dropdown-item" href="#a10.3">10.3 代谢通路图</a>
                </div>
            </div>
"""


pathway_content = """
            <hr>
            <div id="a10" class="anchor"></div>
            <h2>十 代谢通路分析（<a href="09-PathwayTopoEnrichment">09-PathwayTopoEnrichment目录</a>）</h2>

            <div id="a10.1" class="anchor"></div>
            <h3>10.1 富集分析</h3>
            <p>富集分析的目的寻找在某个生物学过程中起关键作用的生物通路, 从而揭示和理解生物学过程的基本分子机制。对于代谢组数据，在富集分析之前，我们需要挑选一些重点关注的代谢物，通常是那些在分组间有显著差异（T检验，p<0.05）的代谢物，再看这些代谢物都在哪些代谢通（KEGG物种特异性代谢通路，不同物种的代谢通路略有区别）路中出现，通过计算这些代谢通路的ORA（Over-Representation Analysis，过表达分析）p值，从而判断关注的代谢物（显著差异的代谢物）是否在这些代谢通路中显著富集。图10-1中显示了差异代谢物显著富集的代谢通路。这些代谢通路可能在所研究的生物学过程中有重要意义。</p>


            <div class="image-block">
                <img src="FiguresTablesForReport/图10-1.png">
                <h6>图10-1 ORA富集分析</h6>
                <p class="image-description">横坐标为富集倍数，它是代谢通路中， 观测代谢物数/理论代谢物数。p值的大小用颜色表示，颜色越深，p值越小。 </p>
            </div>


            <div id="a10.2" class="anchor"></div>
            <h3>10.2 拓扑分析</h3>

            <p>即便关注的代谢物在某个通路中显著富集，我们仍然不知道这些代谢物在这个代谢通路中是否起到关键作用，代谢物对代谢通路究竟有多大影响（Impact）。拓扑分析能够计算关注的代谢物在代谢通路中的作用大小（用Impact衡量）。例如，如果代谢通路中，一个代谢物下游没有其他任何代谢物或基因，那么我们可以判断这个代谢物对代谢通路的影响基本为零。相反，如果一个代谢物非常靠近上游，且下游有很多其他代谢物和基因，那么这个代谢物对代谢通路的影响一定很大。我们通常将拓扑分析和富集分析结合起来，从而判断一个代谢通路是否在所研究的生物学过程中起到关键作用。图10-2中蓝色区域的代谢通路是ORA富集分析中显著的代谢通路，纵坐标展示了这些代谢通路在拓扑分析中的Impact。</p>


            <div class="image-block">
                <img src="FiguresTablesForReport/图10-2.png">
                <h6>图10-2 ORA富集分析和拓扑分析</h6>
                <p class="image-description">横坐标为ORA分析p值，蓝色区域是显著的（p<0.05）；纵坐标为拓扑分析Impact。 </p>
            </div>

            <div id="a10.3" class="anchor"></div>
            <h3>10.3 代谢通路图</h3>
            <p>代谢通路图能够比较直观的反映代谢物的上下游关系，作用模式，以及代谢通路的拓扑结构，同时能找到与代谢物关联的基因。图10-3是只包含代谢物的代谢通路，红色代谢物是关注的代谢物，即在分组间有显著差异的代谢物。图10-4是同时包含代谢物和基因的代谢通路（KEGG完整的代谢通路），有颜色的代谢物是在分组间有显著差异的代谢物，颜色对应分组，表示在对应分组中代谢物含量较高（相对其它分组）。</p>


            <div class="image-block">
                <img src="FiguresTablesForReport/图10-3.png">
                <h6>图10-3 只包含代谢物的代谢通路</h6>
                <p class="image-description">红色代谢物是在分组间有显著差异的代谢物。 </p>
            </div>

            <div class="image-block">
                <img src="FiguresTablesForReport/图10-4.png">
                <h6>图10-4 代谢物和基因的代谢通路</h6>
                <p class="image-description">有颜色的代谢物是在分组间有显著差异的代谢物，颜色对应分组，表示在对应分组中代谢物含量较高（相对其它分组）。 </p>
            </div>
"""


bar_content = """
            <p>在所有代谢物中，有的代谢物是在生物体内扮演着特定角色的，比如激素，维生素等。我们将所有代谢物用<a href="https://www.genome.jp/kegg-bin/get_htext?br08001">KEGG数据库br08001</a>进行注释，得到代谢物所扮演的生物学角色，然后统计每个生物学角色的百分比含量，绘制百分比含量堆积柱形图，如图4-2所示。</p>

            <div class="image-block">
                <img src="FiguresTablesForReport/图4-2.png">
                <h6>图4-2 扮演生物学角色的代谢物百分比堆积柱形图</h6>
                <p class="image-description">横坐标是样品名，根据分组顺序排序，同时用不同颜色标注了不同的分组样本。纵坐标表示各个生物学角色的百分比含量。</p>
            </div>
"""


exp_untarget = """
            <p>该实验中我们利用基于液质联用（LC-MS/MS）的方式研究了样本的代谢组。</p>
            <ol>
                <li>
                    通过Proteowizard软件（v3.0.8789）将获得的原始数据转换成mzXML格式（xcms输入文件格式）
                </li>
                <li>
                    利用R（v3.1.3）的XCMS程序包进行峰识别（peaks identification）、峰过滤（peaks filtration）、峰对齐（peaks alignment）
                </li>
                <li>
                    得到包括质核比（mass to charge ratio, m/z）和保留时间（retention time）及峰面积（intensity）等信息的数据矩阵；进而得到正负离子模式下的前体分子，导出数据至excel进行后续分析。
                </li>
            </ol>
            <div class="image-block">
                <img src="FiguresTablesForReport/src/image/2-1.png">
                <h6>图2-1 代谢组实验流程</h6>
                <p class="image-description"></p>
            </div>
"""

exp_target = """
            <p>样品检测是基于超高效液相色谱串联质谱仪(UPLC-MS/MS)。具体参见方法说明。</p>
"""


def render_html(in_fp, out_fp=False, **kwargs):
    out_fp = out_fp or in_fp
    with open(in_fp, 'r') as f:
        out = f.read()
        for k, v in kwargs.items():
            out = out.replace("{{%s}}" % (k), v)
    with open(out_fp, 'w') as f:
        f.write(out)


if report_type == "target":
    render_html(in_fp=html_path, qc_content="", pathway_content="", pathway_nav="", proc=exp_target,
                bar_content="", od1="十", od2="十一", od3="十二", od4="图3-3", od5="图4-2")
else:
    render_html(in_fp=html_path, qc_content=qc_content, pathway_content=pathway_content, proc=exp_untarget,
                pathway_nav=pathway_nav, bar_content=bar_content, od1="十一", od2="十二", od3="十三", od4="图3-4", od5="图4-3")
