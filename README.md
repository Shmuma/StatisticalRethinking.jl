# StatisticalRethinking


| **Project Status**                                                               |  **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
|![][project-status-img] | [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] [![][codecov-img]][codecov-url] |

## Introduction

This package contains the Julia versions of the snippets contained in the R package "rethinking" associated with the book [Statisticasl Rethinking](https://xcelab.net/rm/statistical-rethinking/) by Richard McElreath.

## Layout of the package

Instead of having all snippets in a single file, the snippets are organized by chapter and grouped in clips by related snippets. E.g. chapter 0 of the R package has snippets 0.1 to 0.5. I have divided those in 2 clips:

1. `clip_00_01_03.jl` - contains snippets 0.1 through 0.3
2. `clip_00_04_05.jl` - contains snippets 0.4 and 0.5.

These 2 files are in chapters/00. These files are later on process using Literate.jl to create 2 derived versions, e.g. from `clip_00_01_03.jl` in chapters/00:

1. `clip_00_01_03.md` - which is stored in docs/src and included in the documentation
2. `clip_00_01_03.ipynb` - stored in the notebooks directory for use in Jupyter

## Acknowledgements

Richard Torker has taken the lead in developing the Turing versions of the models in chapter 8.

The TuringLang team and #turing contributors on Slack have been extremely helpful!

The mcmc components are based on:

1. [TuringLang](https://github.com/TuringLang)
2. [StanJulia](https://github.com/StanJulia)

## Documentation

- [**STABLE**][docs-stable-url] &mdash; **documentation of the most recently tagged version.**
- [**DEVEL**][docs-dev-url] &mdash; *documentation of the in-development version.*

## Questions and issues

Question and contributions are very welcome, as are feature requests and suggestions. Please open an [issue][issues-url] if you encounter any problems or have a question. 

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://stanjulia.github.io/StatisticalRethinking.jl/latest

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://stanjulia.github.io/StatisticalRethinking.jl/stable

[travis-img]: https://travis-ci.org/StanJulia/StatisticalRethinking.jl.svg?branch=master
[travis-url]: https://travis-ci.org/StanJulia/StatisticalRethinking.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/whhifxtx8jb2208f?svg=true
[appveyor-url]: https://ci.appveyor.com/project/StanJulia/StatisticalRethinking-jl

[codecov-img]: https://codecov.io/gh/StanJulia/StatisticalRethinking.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/StanJulia/StatisticalRethinking.jl

[issues-url]: https://github.com/StanJulia/StatisticalRethinking.jl/issues

[project-status-img]: https://img.shields.io/badge/lifecycle-experimental-orange.svg