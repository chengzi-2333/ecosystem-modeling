#import "@preview/easy-paper:0.2.0": *


#show: project.with(
  title: [对种群动态的数学手段模式化分析方法建立暨对实际农作物生产场景预测分析],
  author: array(("3", "45")),
  abstract: [生态学在农业生产中常被用于预测农作物产量。本文使用计算机分析方法，就不同生态因素对叶菜类农作物产量影响进行建模，得到了综合种内、种间竞争、天敌影响、土壤肥力等生态因素综合作用下的叶菜类农作物种群数量模型，最终得以精准预测、精准干预。],
  keywords: ("农作物产量预测", "计算机模拟", "生态学模型"),
  config: (
    text-size: 12pt,
    title-size: 20pt,
  ),
)

#show: rest => columns(2, rest)

农业是我国的基础产业,在我国有悠久的历史. 新中国成立后,我国人口迅速增长,
但受客观因素影响,我国存在人均占有耕地面积少、耕地品质不足等问题,
农业产能难以满足人口增长的需要.
农业生产中,作物的产量往往是人们亟需关注的问题.
在人教版生物选择性必修二《种群的数量特征》一节中,
我们初步学习了"J形增长曲线"与"S形增长曲线".
为进一步系统、深入地学习,并应用于实际生活中,
我们开展了本次研究,结合生态学知识与信息技术,
将实际生产中难以量化的生态学因素,
使用数学语言精确表达,建立数学模型,使用计算机进行数值分析.
考虑农业经济效益影响,我们对农田生态系统进行建模,
并结合现实环境因素进行分析.

= 常见生态学模型
<常见生态学模型>
== Logistic 增长模型
<logistic-增长模型>
皮埃尔-弗朗索瓦·韦吕勒（P.F.Verhulst, 1804-1849）于 1838 年提出
@logistic_origin, 为描述资源有限环境下种群增长的最简单模型,
也是多数生态学模型的基础. 我们熟知的"S形增长曲线"即为该模型.

=== 数学模型描述
<数学模型描述>
由于生境中资源有限,考虑引入环境对种群的影响,
可使用比值$frac(k - x, k)$限制种群增长,
并添加常数$r$以控制增长速率.如此,可以得到
$ #box(stroke: black, inset: 3pt, [$ frac(upright(d) x, upright(d) t) = r x (1 - x / k) . $]) $
其中$x$可表示种群数量或种群密度,$t$为时间, $k$为生物的, $r$为.
内禀增长率$r$是一个经验数值,可由近似公式 $ r = frac(ln R_0, T) $ 给出.
其中$R_0$为净生殖率,即一个世代中种群的增长倍数,$T$为世代时间.

为一阶常微分方程,通解为
$ x (t) = frac(k, 1 + (k / x_0 - 1) e^(- r t)) , $
函数图像如#link(<fig:logistic>)[1];所示,其中 $x_0 = 3$, $r = 1$,
$k = 100$.

#figure(image("./figures/logistic.pdf"), caption: [
  Logistic 增长模型
])
<fig:logistic>

Logistic 方程十分简洁,但未考虑种间影响. 接下来,我们将在 Logistic
增长模型的基础上,添加种间竞争的影响.

== Lotka-Volterra 竞争模型
<sec:models:lv-comp>
二十世纪二十年代,阿尔弗雷德·洛特卡(Alfred Lotka, 1880-1949)
与维多·沃尔泰拉(Vito Volterra, 1860-1940)
分别独立提出了描述种间竞争的简单数学模型,即,
为生态学中描述种间竞争的基本模型.

=== 数学模型描述
<sec:models:lv-comp:desc>
生物竞争的本质是对有限资源的争夺. 考虑种群之间的生境重叠,
可使用$x_j / k_i$描述第$j$种群对第$i$种群生境的影响.

首先考虑两个种群的竞争
$
  frac(upright(d) x_1, upright(d) t) = r_1 x_1 (1 - x_1 / k_1 - alpha_(1 , 2) x_2 / k_1) ,\
  frac(upright(d) x_2, upright(d) t) = r_2 x_2 (1 - x_2 / k_2 - alpha_(2 , 1) x_1 / k_2) .
$

与相比,
仅添加$alpha_(i , j) x_j / k_i$一项(简称"交互项"),
其中$alpha_(i , j)$为, 描述.

更一般地,若存在多个种群,则可叠加不同种群间的交互项
$
  #box(stroke: black, inset: 3pt, [$ frac(upright(d) x_i, upright(d) t) = r_i x_i (1 - 1 / k_i sum_j alpha_(i , j) x_j) . $])
$

将 Logistic
中$x_i / k_i$一项并入求和,
既是为了保持方程形式简洁,也为引入 提供方便.

=== 交互矩阵
<sec:models:lv-comp:inter-mat>
如前文所述,交互系数$alpha_(i , j)$描述第$j$种群对第$i$种群的影响.
注意到这里若要指定具体的$alpha_(i , j)$,需要$i$, $j$两个值,
若种群数量为$n$,则交互系数共有$n^2$个.
随种群数量$n$增长,大量交互系数若采用列表表示,将难以维护.

此时可以令$i$,$j$分别表示行索引与列索引,
构造一个$n times n$的二维"数表",即
$
  bold(alpha) := mat(delim: "[", alpha_(1 , 1), alpha_(1 , 2), dots.h.c, alpha_(1 , j), dots.h.c, alpha_(1 , n); alpha_(2 , 1), alpha_(2 , 2), dots.h.c, alpha_(2 , j), dots.h.c, alpha_(2 , n); dots.v, dots.v, dots.down, dots.v, dots.down, dots.v; alpha_(i , 1), alpha_(i , 2), dots.h.c, alpha_(i , j), dots.h.c, alpha_(i , n); dots.v, dots.v, dots.down, dots.v, dots.down, dots.v; alpha_(n , 1), alpha_(n , 2), dots.h.c, alpha_(n , j), dots.h.c, alpha_(n , n)) ,
$
如此,交互系数便可简单地使用单个符号$bold(alpha)$表示.

#figure(image("./figures/lotka-volterra-competition.pdf"), caption: [
  Lotka-Volterra 竞争模型
])
<fig:lv-comp>

#link(<fig:lv-comp>)[2];为
$upright(bold(x_0)) = angle.l 25 , 50 angle.r$,
$upright(bold(r)) = angle.l 1 , 1 angle.r$,
$upright(bold(k)) = angle.l 100 , 150 angle.r$,
$bold(alpha) = mat(delim: "[", 1, 2; 3, 1)$ 时方程的解.

Lotka-Volterra 竞争模型描述了种间竞争关系.
接下来,我们将进一步引入捕食关系.

== Lotka-Volterra 捕食模型
<sec:models:lv-pred>
经常用于描述生物系统中,
掠食者与猎物进行互动时的动态过程,也就是两者族群规模的消长.
同样地,该模型分别由阿弗雷德·洛特卡与维多·沃尔泰拉于 1925 年与 1926
年独立发表.

=== 数学模型描述
<数学模型描述-1>
首先考虑一对捕食者与猎物. 考虑捕食者对猎物增长率的影响,
可在猎物的增长率$r_1$上添加一项$- beta_(1 , 2) x_2$,
并设$beta_(1 , 2) > 0$,有
$ frac(upright(d) x_1, upright(d) t) = x_1 (r_1 - beta_(1 , 2) x_2) . $
我们暂时设$r_1 > 0$,这表示我们假设猎物在无捕食者的情况下($beta_(1 , 2) = 0$)呈指数增长.

对于捕食者,要添加$beta_(2 , 1) x_1$一项,因为在食物充足情况下有助于其增长,
于是有
$ frac(upright(d) x_2, upright(d) t) = x_2 (- r_2 + beta_(2 , 1) x_1) . $
注意到此处$r_2$前为负号,这表示捕食者在食物匮乏情况下($beta_(2 , 1) arrow.r 0$),
其种群数量有下降趋势.

观察与的形式,
发现系数前符号可并入常数$r$, $beta$,于是可以写出一般形式
$ frac(upright(d) x_i, upright(d) t) = x_i (r_i + beta_(i , j) x_j) . $
综合多个捕食关系的影响,我们得到
$ #box(stroke: black, inset: 3pt, [$ frac(upright(d) x_i, upright(d) t) = x_i (r_i + sum_j beta_(i , j) x_j) . $]) $
由于的$beta_(i , j)$与#link(<sec:models:lv-comp:desc>)[1.2.1];中的$alpha_(i , j)$一样,
均表示了, 因此我们也将$beta_(i , j)$称为,
且由于形式的一致性,同样可以构建由$beta_(i , j)$构成的
$
  bold(beta) := mat(delim: "[", beta_(1 , 1), beta_(1 , 2), dots.h.c, beta_(1 , j), dots.h.c, beta_(1 , n); beta_(2 , 1), beta_(2 , 2), dots.h.c, beta_(2 , j), dots.h.c, beta_(2 , n); dots.v, dots.v, dots.down, dots.v, dots.down, dots.v; beta_(i , 1), beta_(i , 2), dots.h.c, beta_(i , j), dots.h.c, beta_(i , n); dots.v, dots.v, dots.down, dots.v, dots.down, dots.v; beta_(n , 1), beta_(n , 2), dots.h.c, beta_(n , j), dots.h.c, beta_(n , n)) .
$

特别指出,交互矩阵的对角线元素$beta_(i , i)$应满足$beta_(i , i) lt.eq 0$,
否则右侧将出现平方增长项$beta_(i , i) x_i^2$,导致该种群无限增长.

#figure(image("./figures/lotka-volterra-predation.pdf"), caption: [
  Lotka-Volterra 捕食模型
])
<fig:lv-pred>

#link(<fig:lv-pred>)[3];为
$upright(bold(x_0)) = angle.l 60 , 20 angle.r$,
$upright(bold(r)) = angle.l 0.07 , - 0.1 angle.r$,
$bold(beta) = mat(delim: "[", 0, - 0.001; 0.001, 0)$ 时方程的解.

在引入三种基本模型后,我们已初具构建通用模型的条件.
下一节中,我们将进行一些收尾工作,整合三种基本模型的结论,
最终给出通用模型的表达式.

== 通用模型
<通用模型>
通用模型(#strong[G];eneral #strong[M];odel, #strong[GM];),
也可称为大统一模型(#strong[G];rand unified #strong[M];odel),
由三种基本生态学模型导出,取特定参数时可退化为三种基本模型.
它是研究复杂生态系统演化的基础, 在中发挥重要作用.

在接下来的几节中,我们将基于三种基本模型,引入现实因素,逐步构建通用模型.

=== 种间交互
<种间交互>
在#link(<sec:models:lv-comp:inter-mat>)[1.2.2];中,
我们引入了交互矩阵$bold(alpha)$, $bold(beta)$,
这使得我们可以使用更紧凑的形式表示交互系数. 但交互矩阵的

=== 周期性变化
<周期性变化>
在生态系统中,种群可能随环境演化呈周期性变化,例如季节变化.
为引入周期性变化,需要引入周期函数,如正弦型函数$a sin (omega t + phi.alt)$

=== 随机化
<随机化>
周期性变化相对准确地描述了环境变化.
然而,现实环境中存在大量不确定因素,例如难以预测的极端事件.
确定性模型难以描述随机事件, 因此,我们需要引入随机化函数.

考虑到极端事件与常见事件发生概率的差异, 可以使用常见的
描述随机事件对种群的影响 $ X tilde.op N (mu , sigma^2) . $
其中$X$为$mu$为均值,表示随机变量, $sigma$为标准差,描述随机变量.

=== 最终形式
<最终形式>
基于上述分析,我们可以得到通用模型的最终形式
$
  #box(stroke: black, inset: 3pt, [$ #box(stroke: black, inset: 3pt, [$ frac(upright(d) x_i, upright(d) t) = x_i (1 - sum_j alpha_(i , j) x_j) (r_i + sum_j beta_(i , j) x_j) + f_i (t) , $]) $])
$
或更为紧凑的矩阵形式
$
  #box(stroke: black, inset: 3pt, [$ #box(stroke: black, inset: 3pt, [$ frac(upright(d) upright(bold(x)), upright(d) t) = upright(bold(x)) circle.stroked.tiny (upright(bold(I)) - bold(alpha) upright(bold(x))) circle.stroked.tiny (upright(bold(r)) + bold(beta) upright(bold(x))) + upright(bold(f)) . $]) $])
$

#figure(image("./figures/general.pdf"), caption: [
  通用模型
])
<fig:general>

#link(<fig:general>)[4];为
$upright(bold(x_0)) = angle.l 60 , 20 angle.r$,
$upright(bold(r)) = angle.l 0.09 , - 0.1 angle.r$,
$upright(bold(k)) = angle.l 100 , 150 angle.r$,
$bold(alpha) = mat(delim: "[", 1, 2; 3, 1)$,
$bold(beta) = mat(delim: "[", 0, - 0.001; 0.001, 0)$,
$X tilde.op N (0 , 3^2)$ 时方程的解.

= 实例分析
<sec:example>
本章节中,我们将基于通用模型,结合两个具体实例进行分析.

== 稻蟹共生
<sec:example:rice-crab>
稻蟹共生是一种将水稻种植与河蟹养殖相结合的生态农业模式,
通过生物间的互利关系,实现生态效益与经济效益的双赢.
在此生态系统中,我们考虑了杂草、水稻、田间害虫, 形成如下图所示的食物网

= 结论与展望
<结论与展望>
== 数学与农田的对话
<数学与农田的对话>
数学建模作为现实生活中的重要工具, 使用简洁的语言描述事物及其演化.
在作物产量的任务中, 尤其是从理论模型到实际应用的场景,
我们仍然面临诸多挑战.

== 计算机模拟与实际应用的作用
<计算机模拟与实际应用的作用>
=== 预见性
<预见性>
原始的农业生产在有限的经验下,
需要克服不确定的气候、资源等因素,才能提升产量.
计算机的模拟可以使经验判断成为定量计算.
值得注意的是,这种预见性并不要求模型的绝对准确.
复杂计算机模型的校准需要明确处理结构的误差——意识到模型只是对现实的近似,
结果本身包含一定误差,这对决策者而言也是非常重要的信息.

=== 优化性
<优化性>
如上文所述,模型只是对现实的近似. 显而易见,模型越逼近现实,其实用价值越高.
计算机模拟与优化算法的结合,为解决这类问题提供了可行路径.
优化算法可以自我改进深度优化算法,探索许多人类无法察觉的微小因素,
从而做出更加精准的判断.

=== 普及性
<普及性>
理想情况下,每一个农民都要掌握一定的气象、经济学知识,
能准确解读天气数据,能优化资源配置,但现实不是这样的.
计算机模拟在降低专业知识应用门槛方面可以发挥独特作用.
农民不必理解模型的内部细节,不必掌握编程技能就能获得准确可靠的数据.

=== 总结
<总结>
数学建模与计算机模拟的价值,终究要由实践来检验.
在这个意义上,模型永远只是模型——是对现实的简化与近似. 但这种简化如果得当,
就能帮助我们从复杂中识别简单,从噪声中辨别信号,从不确定性中锚定决策.
这或许正是"对种群动态的数学手段模式化分析"对类似稻蟹共生场景的最根本贡献.
当然了,我们也要看到研究的局限性. 从空间的角度看,如在 LV
竞争模型中,如果我们拥有一个更加广阔的充足空间,
就需要考虑到种群的空间分布对于其数量变化的影响.
例如大型食肉动物对小型动物的威慑作用.
从社会的角度看,计算机模拟的能力仍然不足.
例如从人类影响的程度上看,就本研究的例子而言,当我们进一步思考实际意义,
我们需要考虑到人工对害虫的扑杀和人工施肥等因素. 我们相信在不远的将来,
计算机技术可以成熟的被运用到实际生产生活当中,让人类的生活更加美好.

#bibliography("reference.bib")
