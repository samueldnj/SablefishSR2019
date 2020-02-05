
# Renders the SR into body, tandf and appendices

# Body into pdf and docx
bookdown::render_book(  input = "index.Rmd", 
                        clean = TRUE, 
                        output_format = "csasdown::sr_pdf", 
                        config_file = "_bookdown.yml")

bookdown::render_book(  input = "index.Rmd", 
                        clean = TRUE, 
                        output_format = "csasdown::sr_word", 
                        config_file = "_bookdown.yml")

# TandF into pdf and docx
bookdown::render_book(  input = "index.Rmd", 
                        clean = TRUE, 
                        output_format = "csasdown::sr_pdf", 
                        config_file = "_tandf.yml")
bookdown::render_book(  input = "index.Rmd", 
                        clean = TRUE, 
                        output_format = "csasdown::sr_word", 
                        config_file = "_tandf.yml")

# Appendices
bookdown::render_book(  input = "index.Rmd", 
                        clean = TRUE, 
                        output_format = "csasdown::sr_pdf", 
                        config_file = "_appendices.yml")
bookdown::render_book(  input = "index.Rmd", 
                        clean = TRUE, 
                        output_format = "csasdown::sr_word", 
                        config_file = "_appendices.yml")