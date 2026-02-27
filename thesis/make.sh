#!/bin/bash

git clean -fX
latexmk -xelatex article.tex
git clean -fX
