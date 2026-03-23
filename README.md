# 对种群动态的数学手段模式化分析方法建立分析模式暨对实际农作物生产场景预测分析

## $\LaTeX$ 编译

### 编译 `.tex` 文件

```sh
latexmk -xelatex <文件名>
```

### 清理临时文件

```sh
latexmk -bibtex -c  # 小写 "c"
```

### 清理所有 (包括pdf) 编译产物

```sh
latexmk -bibtex -C  # 大写 "C"
```
