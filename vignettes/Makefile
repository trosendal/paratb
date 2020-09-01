model.html: clean
	Rscript -e "library('rmarkdown');render('model.Rmd')"
	cp model.html ../inst/doc

clean:
	rm -rf model_cache
	rm -rf model_files
	rm -f model.html
