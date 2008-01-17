# Declarations
LaTeX = pdflatex      # Latex program
BibTeX = bibtex       # Bibtex program
reflib = mdlrefs
target := manual

#ifndef endnot
#endnot = >/dev/null
#endif

.PHONY: default all clean

default: all

# Defines some variables
PDF := $(foreach targ,$(target),$(targ).pdf)
BBL := $(foreach targ,$(target),$(targ).bbl)
MYT := $(foreach targ,$(target),$(targ).myt)

#all: $(target).pdf
all: 
	make -s $(PDF)

#$(target).myt: $(target).tex
#ifneq ($(shell if [ -e $(target).myt ]; then echo 1; fi),1)
#	touch $(target).myt
#	$(LaTeX) $(target).tex $(endnot)
#endif

#$(target).bbl: $(target).myt mdlrefs.bib
#	$(BibTeX)  $(target) $(endnot)

#$(target).pdf: $(target).tex $(target).bbl
#	$(LaTeX) $(target).tex
#	$(LaTeX) $(target).tex

# Patern rules
# Syntax is
# targets ...: target-pattern: dep-patterns ...
#	commands
#	...

$(MYT): %.myt: %.tex
#ifneq ($(shell if [ -e $*.myt ]; then echo 1; fi),1)
#	touch          $@
#	$(LaTeX) $< $(endnot)
#endif
	if [ ! -e $*.myt ]; then touch $@; $(LaTeX) $< $(endnot); echo "$(LaTeX) $<"; fi

$(BBL): %.bbl: %.myt $(reflib).bib
	echo "$(BibTeX)  $*"
	$(BibTeX)  $* $(endnot)

$(PDF): %.pdf: %.tex %.bbl
#ifneq ($(shell if [ $*.aux -nt $*.bbl ]; then echo 1; fi),1)
#	$(LaTeX) $*.tex $(endnot)
#endif
	if [ ! $*.aux -nt $*.bbl ]; then $(LaTeX) $*.tex $(endnot); echo "$(LaTeX) $<"; fi
	$(LaTeX) $*.tex $(endnot)
	echo "$(LaTeX) $*.tex"

# Implicit rules (suffix rules)
# NB. $@ = current target, $< = implied source, $* = filename without suffix.
# Syntax is
# .insuffix.outsuffix:
#	command
#.tex.myt:
#ifneq ($(shell if [ -e $@ ]; then echo 1; fi),1)
#	touch $@
#	$(LaTeX) $< $(endnot)
#endif
#
#.tex.pdf:
#ifneq ($(shell if [ $*.aux -nt $*.bbl ]; then echo 1; fi),1)
#	$(LaTeX) $< $(endnot)
#endif
#	$(LaTeX) $< $(endnot)
#
#.aux.bbl:
#	$(BibTeX)  $* $(endnot)
#
#.tex.aux:
#	$(LaTeX) $< $(endnot)

#.bbl: $*.myt $(reflib).bib

#.pdf: $*.tex $*.bbl

#bib_target: $(target).bbl
#	$(BibTeX) $(target)
#
#tex_target: $(target).aux
#	$(LaTeX) $(target).tex
##	$(BibTeX) $(target)
##	$(LaTeX) $(target).tex
##	$(LaTeX) $(target).tex
#
#all_tex: tex_target bib_target 
#	$(LaTeX) $(target).tex

clean:
	rm -vf *.aux *.bbl *.blg *.log *.toc
	#rm -vf $(target).dvi $(target).ps $(target).pdf $(target).myt
	rm -vf $(foreach targ,$(target),$(targ).dvi) $(foreach targ,$(target),$(targ).myt)
	rm -vf $(foreach targ,$(target),$(targ).ps) $(foreach targ,$(target),$(targ).pdf)
	rm -vf *~

# For the implicit rules to work we must declare the suffixes used here
.SUFFIXES: .tex .aux .toc .bbl .blg .log .dvi .ps .pdf .myt
