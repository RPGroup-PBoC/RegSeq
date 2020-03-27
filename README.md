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

#### **Installing module**
In order to use the functions within the `RegSeq` module it is necessary to
install the package locally. This can simply be done by navigating in the
terminal to the main project directory and typing the command
```
pip install -e ./
```
The 'setup.py' file will take care of the installation. The `-e` option within
the package allows for the constant update of the package as it might be
subject to changes as the project is being developed.

The modules contained in the package include:

`module1.py` : explanation here.

`module1.py` : explanation here.

# License Information
<img src="https://licensebuttons.net/l/by-nd/3.0/88x31.png"> This work is
licensed under [CC-BY-ND](https://creativecommons.org/licenses/by-nd/4.0/). All
software is issued under the standard MIT license which is as follows:

```
Copyright 2019, The authors

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
