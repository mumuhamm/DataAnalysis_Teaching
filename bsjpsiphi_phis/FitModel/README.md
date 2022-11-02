
# FitModel
Fitter for the JPsi/Phi angular analysis.
Developers: M.D. Alibordi, A. Boletti, A. Bragagnolo, G. Fedi, E. Lusiani.

## Files in  this folder
### Main Fit Macros
* `Pwave_Reco_Model_recotag.C`
    Launch macro for the fitting reco samples, calls all the single fit routines
 * `Pwave_Gen_Model_closure.C`
  Launch macro for fitting a gen samples
* `Pwave_Model.C`
  Main fitting code, should be called from  `Pwave_Reco_Model_recotag.C` or `Pwave_Gen_Model_closure.C`
  
### Single Fit Routines
These macros can be run standalone and their output can be use later in the full fit. These macro will save the results inside a file in a directory called `ouputs`.

* `Bs_angEff_proj.C(TString recofile, TString genFile)`
    Macro to compute the angular function using Legendre polynomials. Also computes the angular efficiency which is used inside the function. From Alessio Boletti work.
* `Bs_ctEff_fit.C(TString recofile, TString genFile)`
    Macro to compute and fit the ct efficiency
* `Bs_ctErrPdf_fit.C(TString recofile)`
    Macro to fit the ct error pdf
* `Bs_massPdf_fit.C(TString recofile)`
    Macro to fit the mass pdf
    
### Closure Tests
These macros perform validation tests.
* `Bs_angEff_closure.C`
  Macro to check if the angular efficiency functions describe well the angular variables.

### Utilities
* `utils.h`
    Common utils to all the macro, like getting a tree from a file or saving a canvas in multiple formats.
* `definitions.h`
    Definitions of the observables used in the fit, needed to unify them over all the macros
* `Pwave_Model_plot.C`
  Plotting code


## How to run the gen fit
Simply call the macro with the gen sample as argument

```root -l -b -q 'Pwave_Gen_Model_closure.C("path_to_gen/ntuBsGEN.root", false)'```

The second argument is a boolean flag to switch plotting on/off (takes a lot of time currently). Default is ON.

## How to run the reco fit
The reco fit macro can be used in two ways:
### Standalone execution

 1. Only reco/gen datasets path and the initial parameter flag needed, all subroutines will be called automatically
```
root -l -b -q 'Pwave_Reco_Model_recotag.C("path_to_reco/fittree_ntuBsMC2017.root", "path_to_gen/ntuBsGEN.root", isNotDG0, "fit_suffix")'
```
### Via subroutines
This allows to runs only certain parts of the fit or to avoid running parts that do not need to be refitted (e.g. using eff evaluated in sampleA to fit sampleB).

 1.  Run one of more subroutine to produce the needed fit/efficiencies files. This will produces files in a directory called `outputs`

```
root -l -b -q 'Bs_ctErrPdf_fit.C("path_to_reco/fittree_ntuBsMC2017.root")'
root -l -b -q 'Bs_massPdf_fit.C("path_to_reco/fittree_ntuBsMC2017.root")'
root -l -b -q 'Bs_angEff_proj.C("path_to_reco/fittree_ntuBsMC2017.root", "path_to_gen/ntuBsGEN.root")'
root -l -b -q 'Bs_ctEff_fit.C("path_to_reco/fittree_ntuBsMC2017.root", "path_to_gen/ntuBsGEN.root")'
```
   
 2.  Run the main fit macro, giving the previously produced files as arguments
```
root -l -b -q 'Pwave_Reco_Model_recotag.C("path_to_reco/fittree_ntuBsMC2017.root", "path_to_gen/ntuBsGEN.root", isNotDG0, "outputs/ctEffDir.root", "outputs/angEffDir.root", "outputs/massPdfDir.root", "outputs/ctErrPdfDir.root", false)'
```
last parameter is again the plotting flag
### MultiSample execution
 1. Optional: run the subroutines. Subroutines except for mass have a last parameter which is the suffix which is added to the files/plots before saving. 
 
    Mass fit is always shared between datasets and is best left to the main fit
    
 2. Create a config file containing, in order
    * A label for the dataset
    * The path to the dataset file
    * The path of the ct efficiency file (or null if not available)
    * The path of the ang efficiency file, or null
    * The path of the ct error file, or null
    
    for each dataset in a line
    
    Example:
    
```
y17 path_to_reco/fittree_ntuBsDG0MC2017.root outputs/ctEffDiry17.root null null
y18 path_to_reco/fittree_ntuBsDG0MC2017.root null outputs/angEffDiry18.root outputs/ctErrpdfDiry18.root
```
 3. Run the main fit macro, with the config as argument
```
root -l -b -q 'Pwave_Reco_Model_recotag.C("path_to_config", "path_to_gen/ntuBsGEN.root", isNotDG0, "outputs/massPdfDir.root", wantPlot)'
```
### Toy MCs
1. Run the mistag/tag fit. Mistag and tag are fitted separately, but are needed during generation

    `root -l -b -q 'Bs_tagMistag_fit("path_to_data", "data_label")'`

    where data label is an arbitrary label that will be appended to the produced outputs file

2. Run the generator macro

    Generator needs as inputs:

    * the model file (produced during the main fit)
    * a config file that links the mistag output to the dataset label used in the main fit
    * the path to the (preexisting) output dir where to save the samples
    * the generated sample size
    * the number of samples to generate in this execution (initialization is not too fast and generating more samples with one execution is useful)
    * an integer to be used as seed (same seed = same samples)

    The format of the mistag config file is similar to the one used in the main fit, meaning it is a label followed by the output file name of the mistag/tag fit macro

    Labels must be the same as the ones used in the main fit

    `root -l -q -b 'Pwave_toyMC_generator("path_to_model", "mistag_config_file", "output_dir", sample_size, samples_per_generation, seed)'`

    This routine should be called many times with different seeds, same amount of samples and output dirs in a pattern similar to "prefix_n", where n is an index starting at 0, as that is the format the aggragator expects. 
    The number of samples per generation shouldn't be too high (as they are generated sequentially and occupy RAM), but not too low either (as you would initialize the generator often).
    5/10 is a good compromise.

3. Run the fitter.
    Fitter takes as input

    * the model file from the main fit the dir where the samples are stored
    * the index of the sample inside the dir (this expect tha sample to follow the convention of the generator, which places n files with names data.i.root inside the output directory)
    * the number of cpus

    `root -l -q -b 'Pwave_toyMC_fitter("path_to_model", "generation_dir", sample_index, nCPU)'`

4. Run the aggregator, that collects all the fit files and creates a tree with the results

    `root -l -b -q 'Pwave_toyMC_aggregator("prefix", total_number_of_generations, samples_per_generation, "path_to_model")'`

5. Run the visualizer, to produce the plots, passing in the name of the output file from the aggregator.

    `root -l -b -q 'Pwave_toyMC_visualizer("outputs/toyMCPars.root")'`