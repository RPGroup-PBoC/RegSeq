import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="regseq",
    version="0.0.4",
    author="Bill Ireland, Tom Roeschinger, Niko McCarty, Rob Phillips",
    author_email="troeschi@caltech.edu",
    description="This repository contains all active research materials for the Reg-Seq project.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/RPGroup-PBoC/RegSeq",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)