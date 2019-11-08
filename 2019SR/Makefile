export PATH := /Library/TeX/texbin:$(PATH)    # add LaTeX path
export PATH := /usr/local/bin:~/.local/bin:$(PATH) # add pandoc-citeproc-preamble

# Cluster targets
pdf: index.pdf tandf.pdf appendices.pdf
docx: index.docx tandf.docx appendices.docx

all: pdf docx

## pdf
index.pdf: index.Rmd Makefile
	Rscript -e 'bookdown::render_book(input = "index.Rmd", clean = TRUE, output_format = "csasdown::sr_pdf", config_file = "_bookdown.yml")'

# docx
index.docx: index.Rmd Makefile
	Rscript -e 'bookdown::render_book(input = "index.Rmd", clean = TRUE, output_format = "csasdown::sr_word", config_file = "_bookdown.yml")'

tandf.pdf: index.Rmd Makefile
	Rscript -e 'bookdown::render_book(input = "index.Rmd", clean = TRUE, output_format = "csasdown::sr_pdf", config_file = "_tandf.yml")'

# docx
tandf.docx: index.Rmd Makefile
	Rscript -e 'bookdown::render_book(input = "index.Rmd", clean = TRUE, output_format = "csasdown::sr_word", config_file = "_tandf.yml")'


appendices.pdf: index.Rmd Makefile
	Rscript -e 'bookdown::render_book(input = "index.Rmd", clean = TRUE, output_format = "csasdown::sr_pdf", config_file = "_appendices.yml")'

# docx
appendices.docx: index.Rmd Makefile
	Rscript -e 'bookdown::render_book(input = "index.Rmd", clean = TRUE, output_format = "csasdown::sr_word", config_file = "_appendices.yml")'

