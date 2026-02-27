#!/bin/bash

latexmk -C
latexmk -xelatex
latexmk -c -bibtex
