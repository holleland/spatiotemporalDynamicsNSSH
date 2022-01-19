Body growth, somatic condition and gonad development of Norwegian spring
spawning herring (Clupea harengus) linked to large-scale spatiotemporal
dynamics in distribution
================

*Contributors: Eydna Í Homrum<sup>1\*</sup>, Guðmundur J.
Óskarsson<sup>2</sup>, Kotaro Ono<sup>3</sup>, Sondre
Hølleland<sup>3</sup>, Aril Slotte<sup>3</sup>*

<sup>1</sup> Faroe Marine Research Institute, Tórshavn, Faroe Islands,
<br> <sup>2</sup> Marine and Freshwater Institute, Hafnafjörður, Iceland
<br> <sup>3</sup> Institute of Marine Research, Bergen, Norway.<br>
<sup>\*</sup> Corresponding author: <eydnap@hav.fo><br>

[Link to published paper goes
here](https://www.frontiersin.org/journals/marine-science)

The purpose of this github repository is to ensure transparency of the
results presented in **Homrum et. al. (2022)**. Here you will find the
necessary R- and C++ code to fit the models and produce most of the
figures in the paper. Due to restrictions detailed below the data is not
included, but may be shared upon request.

## Code

The code is available in the *R/* folder and divided into different
scripts. Each script will save the results in the output folder and this
is then loaded (if needed) in the other scripts. Each script start with
a clean R environment. No result files are uploaded to github, but
should be reproducible by executing the code in the folders (following
the numbering). In short, the *R/0\_data.R* file organizes the data,
*R/1\_analysis.R* fits all the models and stores them in the output
folder. The remaining scripts creates plots based on the data and the
fitted models. The script *R/residual\_diagnostics.R* contains helping
functions for residual analysis.

The models are implemented using *Template Model Builder* (TMB) and the
C++ code can be found in the *src/* folder.

## Data availability

The data used in the **Homrum et. al. (2022)** are stored in a joint
database (PGNAPES), and can only be accessed by the consent of all
parties. Requests to access these datasets should be directed to
Guðmundur J. Óskarsson, <gjos@hafogvatn.is>; Aril Slotte,
<aril.slotte@hi.no>; and Jan Arge Jacobsen,
[janarge@hav.fo](emailto:janarge@hav.fo).

## Author’s github accounts

Kotaro Ono - [kotkot](https://github.com/kotkot)

Sondre Hølleland - [holleland](https://github.com/holleland)

Eydna Í Homrum - [eydnaihomrum](https://github.com/eydnaihomrum)

## Licence

This project is licensed under the GNU GPLv3 License - see
[LICENSE](LICENSE) for details.

## Funding

This study was partly funded by the project “Sustainable multi-species
harvest from the Norwegian Sea and adjacent ecosystems” the Research
Council of Norway (project number 299554).

## References

Homrum, E.Í., Óskarsson, G.J., Ono, K., Hølleland, S., and Slotte, A.
(2022). Body growth, somatic condition and gonad development of
Norwegian spring spawning herring (Clupea harengus) linked to
large-scale spatiotemporal dynamics in distribution. Manuscript under
review in *Frontiers of Marine Science*.
