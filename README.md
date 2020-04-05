# RegSeq

## Overview
This is a GitHub repository for Reg-Seq, an experimental protocol that determines the binding energy matrix for regulatory regions of the genome. The website for the original Reg-Seq paper can be found [here](https://www.rpgroup.caltech.edu/RNAseq_SortSeq/).

**Check out the Wiki tab to see the full experimental protocol for Reg-Seq**.

## Layout
The repository is split into four main directories, many of which have
subdirectories. This structure has been designed to be easily navigable by
humans and computers alike, allowing for rapid location of specific files and
instructions. Within each directory is a `README.md` file which summarizes the
purpose of that directory as well as some examples where necessary. 

`code` contains Jupyter notebooks to perform computational steps of the Reg-Seq protocol from start to finish. Where all of the executed code lives. This includes pipelines, scripts, and
figure files.

`data` contains prior data files, such as designed sequences, negative controls, or output analysis files.

`protocol` contains MarkDown files that outline each step to perform Reg-Seq in total.

`RegSeq` contains Python files which can be easily executed to perform your own analyses.

#### **Package Dependencies**
Some of the code files in this package rely on `mpathic`, a software package for quantitative modeling of massively parallel experiments and developed by [GitHub user jbkinney](https://github.com/jbkinney). A link to the GitHub page for the mpathic package is [available here](https://github.com/jbkinney/mpathic).

Currently, installation of this package requires Linux or Mac OS. To install the mpathic package, please run the following commands within an Anaconda terminal:

`conda create -n mpathic_env python=3.6.9 pip` --> create a Python environment for version 3.6.9

`source activate mpathic_env` --> activate the python 3.6.9 environment

`conda install mpathic -c wireland` --> install the mpathic package, using username 'wireland'

`conda install mpmath` --> install another package, called mpmath

`mpathic learn_model --help` --> verify that the installation proceeded as expected. Running this command should populate the command terminal with a list of available functions.

If you encounter any issues with installation, please contact us through GitHub.

The modules contained in the RegSeq package include:

`create_key_to_match_sequence_to_barcode.py` : initial analysis script to match barcodes to sequences.

`do_many.py` : A script to parallelize analysis via Amazon AWS.

`learn_model_mut.py` :  Infers the effect of mutation on gene expression for each base pair using MCMC. This is crucial for converting to an information footprint.

`peak_utils.py` : functions for information footprints.

`plot_informationfootprint.py` : Python script to generate information footprints from sequencing data.

`utils.py` : functions for analysis.

# Setting up a kernel for MPAthic
After installing the MPAthic package, those notebooks and code blocks with dependencies can be executed. Note that modern versions of Anaconda **do not** automatically create kernels for new environments. Therefore, you must manually create a kernel while in the mpathic_env. To do that, execute:

`python -m ipykernel install --user --name myenv --display-name "Python (myenv)"` 

from the new environment. Then, after opening the Jupyter notebook (e.g. `checkmutations.ipynb`), click on 'Kernel', 'Change Kernel', and select the newly-created kernel. You will now be able to import the package and execute the code.

# License Information
<img src="https://licensebuttons.net/l/by-nd/3.0/88x31.png"> This work is
licensed under [CC-BY-ND](https://creativecommons.org/licenses/by-nd/4.0/). All
software is issued under the standard MIT license which is as follows:

```
Copyright 2020, The authors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
