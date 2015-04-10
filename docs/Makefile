
%.html : %.Rmd
	R --slave -e 'rmarkdown::render("$(firstword $^)")'

%.pdf : %.Rmd
	R --slave -e 'rmarkdown::render("$(firstword $^)")'

analysis.html : analysis.Rmd
