ECRad is a code for the predection of electron cyclotron emission measurements of thermal and non-thermal magnetic confinement fusion plasmas.

## References
At the moment [this paper](https://doi.org/10.1016/j.cpc.2020.107175) is the reference for ECRad.
However, there are two other papers that should be cited alongside the above paper:
* If you use the *primary* model for the absorption coefficent cite [this](https://doi.org/10.1088/0741-3335/49/1/002).
* If you use the *secondary* model for the absorption coefficent cite [this](https://doi.org/10.13182/FST07-A1494).
As per the author's request this reference differs from the reference in the ECRad paper.

You can find the corresponding .bib files in the ``references`` folder. There is no guarantee for the correctness of these refenrence.


## Installation

ECRad_core is released as a conda package from the [ECRad](https://github.com/AreWeDreaming/ECRad) repository, which also holds the conda recipes and the instructions for a local development install.

## Dependencies
Either gfortran or intel ifort and GNU make.

## Usage
ECRad should be exclusively run through its designated (GUI)[https://github.com/AreWeDreaming/ECRad].

## Contributing
At the moment all contributions should be discussed. Pull requests will be welcome once version 1.0 is stable.

## License
[MIT](https://choosealicense.com/licenses/mit/)

Exceptions:

ECRad contains modified versions of fitpack and odepack from the scipy python package which falls under the [BSD license](https://www.scipy.org/scipylib/license.html). 
The corresponding license can be found in license folder.

The source file quadratue.f90 is included in its unmodified form. 
It is published under the [GNU LGPL license](https://choosealicense.com/licenses/lgpl-3.0/mit/). The text of the license can be found in the license folder.
