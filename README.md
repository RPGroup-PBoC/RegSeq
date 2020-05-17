# RegSeq

## Overview
This is a GitHub repository for Reg-Seq, an experimental protocol for deciphering the regulatory architecture of bacterial promoters. Specifically, the goal of Reg-Seq is to take previously uncharacterized or partially characterized promoters and to determine
the constellation of RNAP and transcription factor binding sites and to determine what transcription factors
bind those sites. In addition, the method allows for a determination of energy matrices that characterize,
in k_BT units, the binding strength of transcription factors to their target sites.

The website for the original Reg-Seq paper can be found [here](https://www.rpgroup.caltech.edu/RNAseq_SortSeq/).

**Check out the [Wiki tab](https://github.com/RPGroup-PBoC/RegSeq/wiki) to see the full experimental protocol for Reg-Seq**.

## Installation
To reproduce this work, you will need to use the RegSeq module -- a homegrown Python software package written explicitly for this work. We recommend having installed Anaconda3. The package relies on the `mpathic` package, a software package for quantitative modeling of massively parallel experiments and developed by [GitHub user jbkinney](https://github.com/jbkinney). A link to the GitHub page for the mpathic package is [available here](https://github.com/jbkinney/mpathic). `mpathic` relies on Python 3.6.9. Since this is not the most recent python version, it is very likely that a new Python environment is necessary. To create this environment, simply navigate into the `RegSeq/` folder and run the following line in the terminal:

`conda env create -f environment.yml`

Now there should be new python environment titled `mpathic_env`. The environment can be activated by:

`conda activate mpathic_env`

And deactivated by:

`conda deactivate`

If any commands are run from the command line, it needs to be done inside the created environment. To confirm that the environment is functional, run the following line after activating the new environment:

`mpathic learn_model --help`

Which verifies that the installation proceeded as expected. Running this command should populate the command terminal with a list of available functions.

Finally, to use this environment in a Jupyter notebook, a kernel needs to be created. Therefore, activate the environment and run the following line:

`python -m ipykernel install --user --name mpathic_env --display-name "Python (mpathic)"`

When opening a notebook from the `RegSeq/code` folder, click on 'Kernel', 'Change Kernel', and select the newly-created kernel. You will now be able to import the package and execute the code.

You should be set up to use all code provided in this repository. If you encounter any issues with the installation, please contact us through GitHub.

## Layout
The repository is split into four main directories, many of which have
subdirectories. This structure has been designed to be easily navigable by
humans and computers alike, allowing for rapid location of specific files and
instructions. Within each directory is a `README.md` file which summarizes the
purpose of that directory as well as some examples where necessary.

### notebooks 
Contains Jupyter notebooks to perform computational steps of the Reg-Seq protocol from start to finish. Where all of the executed code lives. This includes pipelines, scripts, and
figure files.

### data 
Contains prior data files, such as designed sequences, negative controls, or output analysis files.


### reqseq
Python module which can be easily executed to perform your own analyses.

### figures
Contains all generated figures.

# License Information
<img src="https://licensebuttons.net/l/by-nd/3.0/88x31.png"> This work is
licensed under a [Creative Commons CC-BY 4.0 Attribution license](https://creativecommons.org/licenses/by/4.0/). All
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
