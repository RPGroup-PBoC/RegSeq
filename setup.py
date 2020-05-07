import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="regseq",
    version="0.0.1",
    author="Bill Ireland, Niko McCarty, Tom Roeschinger, Rob Phillips",
    author_email="nmccarty {at} caltech {dot} edu",
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