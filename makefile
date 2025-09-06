# Usually, only these lines need changing
# No spaces after names
QPAPFILE = SCSK
MAIN = mainMC#mainParallelOptPCs
# QSLIFILE = JMP-update
RDIR = ./Simulations/R
MLDIR = ./Simulations/Matlab
THDIR = ./Theta
QPAPDIR = ./Document
QREPDIR = ./Theta/reports
# QREPDIRSECDIR = ./Theta/reports/sections
# QSLIDIR = ./Quarto-Slides

# Define rho values for Monte Carlo simulations
# RHO_VALUES = 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
RHO_VALUES = 0.0 0.4 0.9 0.1 1.0 0.6 0.3 0.7 0.5 0.8
# list all R files
RFILES := $(wildcard $(RDIR)/*.R)
EXCLUDE := $(wildcard $(RDIR)/_*.R)
# QFILES := $(wildcard $(QPAPDIR)/*.qmd)
QREPFILES := $(wildcard $(QREPDIR)/*.qmd)
# QSECFILES := $(wildcard $(QREPDIRSECDIR)/*.qmd)
# QEXCLUDE := $(wildcard $(QPAPDIR)/_*.qmd)
# RDEP := $(wildcard $(RDIR)/*beta_diff.R)

# excluding files that start with "_" during development
RFILES := $(filter-out $(EXCLUDE),$(RFILES))
# QFILES := $(filter-out $(QEXCLUDE),$(QPAPFILES))
# QFILES := $(filter-out $(QPAPFILE).qmd,$(QPAPFILES))
# RIND := $(filter-out $(RDEP),$(RFILES))

# list all matlab m files
MTLB_FILES := $(wildcard $(THDIR)/*.m)
MTLB_EXCLUDE := $(wildcard $(THDIR)/_*.m)
MTLB_FILES := $(filter-out $(MTLB_EXCLUDE),$(MTLB_FILES))

# Indicator files to show R file has run
OUT_FILES := $(RFILES:.R=.Rout)
# RDEPOUT := $(RIND:.R=.Rout)
OUT_MTLB := $(MTLB_FILES:.m=.log)
OUT_REPS := $(QREPFILES:.qmd=.tex)
MLOUT := $(THDIR)/$(MAIN).log

# Targets


## Default target
main: $(RDIR)/main.Rout #$(filter-out $(RDIR)/main.Rout, $(OUT_FILES))

## Make all
all: matlab paper mail

## Run MATLAB with different rho values
matlab-rho: $(addprefix $(THDIR)/mainMC-r, $(addsuffix .log, $(RHO_VALUES)))

## Run MATLAB with a specific rho value (usage: make matlab-rho-single RHO=0.3)
matlab-rho-single:
	@if [ -z "$(RHO)" ]; then echo "Please specify RHO value: make matlab-rho-single RHO=0.3"; exit 1; fi
	matlab -sd $(THDIR) -logfile $(THDIR)/mainMC-r$(RHO).log -batch 'rho=$(RHO); mainMC'

## Run R files
R: $(OUT_FILES)


mtlb: $(MLOUT) #$(MLDIR)/$(MATLAB).log

theta: $(OUT_MTLB)

report: $(OUT_REPS)

echo:
	echo $(MLOUT) $(OUT_MTLB) $(MTLB_FILES) $(OUT_FILES) $(OUT_REPS)

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

theta_report: $(THDIR)/reports/Theta-report.qmd 
	quarto render $<



## Send mail with results and paper
mail: $(RDIR)/_send_mail.Rout

# Rules
$(RDIR)/%.Rout: $(RDIR)/%.R 
	R CMD BATCH --no-save --no-restore-data $< $@

$(RDIR)/_send_mail.Rout: $(RDIR)/_send_mail.R
	R CMD BATCH --no-save --no-restore-data $< $@

# $(MLDIR)/$(MATLAB).log: $(MLDIR)/$(MATLAB).m
# 	matlab -sd $(MLDIR) -logfile $@ -batch $(notdir $(basename $<))

# $(THDIR)/$(MATLAB).log: $(THDIR)/$(MATLAB).m
# 	matlab -sd $(THDIR) -logfile $@ -batch $(notdir $(basename $<))


$(THDIR)/%.log: $(THDIR)/%.m
	matlab -sd $(THDIR) -logfile $@ -batch $(notdir $(basename $<))

# Rule for running mainMC with specific rho values
$(THDIR)/mainMC-r%.log: $(THDIR)/mainMC.m
	matlab -sd $(THDIR) -logfile $@ -batch 'rho=$*; mainMC'
	
# # Compile main tex file and show errors
$(QPAPDIR)/$(QPAPFILE).pdf: $(QPAPDIR)/$(QPAPFILE).qmd #$(QFILES) #$(OUT_FILES) #$(CROP_FILES)
	quarto preview $<

$(QREPDIR)/%.tex: $(QREPDIR)/%.qmd
	quarto render $<

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
	rm -fv $(THDIR)/mainMC-r*.log

clean-latex:
	rm -fv */*.aux  */*.fdb_latexmk */*.fls */*.lot */*.lof */*.synctex.gz
	rm -fv *.aux  *.fdb_latexmk *.fls *.lot *.lof *.synctex.gz
	rm -fv Document/*/*.aux  Document/*/*.fdb_latexmk Document/*/*.fls Document/*/*.lot Document/*/*.lof Document/*/*.synctex.gz
	rm -fv Theta/*/*.aux  Theta/*/*.fdb_latexmk Theta/*/*.fls Theta/*/*.lot Theta/*/*.lof Theta/*/*.synctex.gz


.PHONY: all clean paper mail mtlb report clean-out clean-ml clean-latex theta theta_report