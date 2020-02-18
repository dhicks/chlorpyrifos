paper.pdf: paper.md
	pandoc paper.md -o paper.pdf -F pandoc-citeproc --pdf-engine=lualatex




