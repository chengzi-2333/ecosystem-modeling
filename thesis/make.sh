#!/bin/bash
git clean -fX
xelatex article
biber article
xelatex article
xelatex article
git clean -fX
