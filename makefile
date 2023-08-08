# Usually, only these lines need changing
QPAPFILE = SCPC
MATLAB = main
# QSLIFILE = JMP-update
RDIR = ./Simulations/R
MLDIR = ./Simulations/Matlab
QPAPDIR = ./Document
# QSLIDIR = ./Quarto-Slides

# list all R files
RFILES := $(wildcard $(RDIR)/*.R)
EXCLUDE := $(wildcard $(RDIR)/_*.R)
# QFILES := $(wildcard $(QPAPDIR)/*.qmd)
# QEXCLUDE := $(wildcard $(QPAPDIR)/_*.qmd)
# RDEP := $(wildcard $(RDIR)/*beta_diff.R)

# excluding files that start with "_" during development
RFILES := $(filter-out $(EXCLUDE),$(RFILES))
# QFILES := $(filter-out $(QEXCLUDE),$(QPAPFILES))
# QFILES := $(filter-out $(QPAPFILE).qmd,$(QPAPFILES))
# RIND := $(filter-out $(RDEP),$(RFILES))


# Indicator files to show R file has run
OUT_FILES := $(RFILES:.R=.Rout)
# RDEPOUT := $(RIND:.R=.Rout)

# MLOUT := $(MLDIR)/$(MATLAB).log
# Targets


## Default target
main: $(RDIR)/main.Rout #$(filter-out $(RDIR)/main.Rout, $(OUT_FILES))

## Make all
all: matlab paper mail

## Run R files
R: $(OUT_FILES)

matlab: $(MLDIR)/$(MATLAB).log

## Make paper
paper: #$(QPAPDIR)/$(QPAPFILE).pdf
	quarto render $(QPAPDIR)/$(QPAPFILE).qmd 
#	quarto render $(QPAPDIR)/$(QPAPFILE).qmd --to pdf -M include-in-header:packages.tex
#	open -a Preview $(QPAPDIR)/$(QPAPFILE).pdf
# pdf:
# 	quarto render $(QPAPDIR)/$(QPAPFILE).qmd --to pdf -M include-in-header:packages.tex
## Make slides
# slides: #$(QSlIFILE).html
# 	quarto render $(QSLIDIR)/$(QSLIFILE).qmd
#	open -a Safari $(QSLIDIR)/$(QSLIFILE).html

## Send mail with results and paper
mail: $(RDIR)/_send_mail.Rout

# Rules
$(RDIR)/%.Rout: $(RDIR)/%.R 
	R CMD BATCH --no-save --no-restore-data $< $@

$(RDIR)/_send_mail.Rout: $(RDIR)/_send_mail.R
	R CMD BATCH --no-save --no-restore-data $< $@

$(MLDIR)/$(MATLAB).log: $(MLDIR)/$(MATLAB).m
	matlab -sd $(MLDIR) -logfile $@ -batch $(notdir $(basename $<))

# # Compile main tex file and show errors
$(QPAPDIR)/$(QPAPFILE).pdf: $(QPAPDIR)/$(QPAPFILE).qmd #$(QFILES) #$(OUT_FILES) #$(CROP_FILES)
	quarto preview $<

# # Compile main tex file and show errors
# $(QSLIFILE).html: $(QSLIFILE).qmd $(OUT_FILES) #$(CROP_FILES)
#     quarto preview "$(QSLIDIR)/$(QSLIFILE).qmd"

# Dependencies
# May need to add something here if some R files depend on others.
$(RDIR)/main.Rout: $(RFILES)

$(RDIR)/100_functions.Rout: $(RDIR)/000_global_vars.R

# $(RDIR)/30_beta_diff.Rout $(RDIR)/40_beta_diff_yoy.Rout $(RDIR)/45_beta_diff.Rout :$(RDIR)/20_functions.R $(RDIR)/15_global_vars.R $(RDIR)/00_reading_data.R

# $(RDIR)/60_plot_discontinuity.Rout $(RDIR)/50_beta_diff.Rout $(RDIR)/30_beta_diff.Rout:$(RDIR)/20_functions.R $(RDIR)/15_global_vars.R $(RDIR)/00_reading_data.R

# $(RDEPOUT) : $(RIND)


# Clean up stray files
clean:
	rm -fv $(OUT_FILES) 
	rm -fv *.Rout *.RData
	rm -fv *.aux *.log *.toc *.blg *.bbl *.synctex.gz *.out *.bcf *blx.bib *.run.xml
	rm -fv *.fdb_latexmk *.fls
#	rm -fv $(TEXFILE).pdf
clean-out:
	rm -fv $(OUT_FILES) *.Rout

clean-ml:
	rm -fv $(MLDIR)/$(MATLAB).log

.PHONY: all clean paper mail matlab