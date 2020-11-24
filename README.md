# PDB-to-GNM
Apply Normal Mode Analysis to a Protein Data Bank structure (PDB) using Gaussian Network Modelling



Last update: 2020-11-19

This repository is administered by [@aarondomenzain](https://github.com/aarondomenzain), Departamento de Física, Área de Física de Líquidos, Universidad Autónoma Metropolitana Unidad Iztapalapa.

## Installation
### Required Software
* [Python 3.x.x] (https://www.python.org/downloads/)
* [PIP](https://pip.pypa.io/en/stable/installing/) 
### Dependencies - Recommended Software
#### Python Libraries
* Numpy
* Matplotlib
* Scipy
* Prody
### Installation Instructions
* Clone `main` branch from [here](https://github.com/aarondomenzain/PDB-to-GNM).
* Open a terminal window and run the following commands:
```
pip install numpy scipy matplotlib prody
```
## Development Guidelines

Anybody is welcome to contribute to the development of this repository, but please abide by the following guidelines.

Each function should start with a commented section describing the function and explaining the parameters. Existing functions can clarify what style should be used. When making *any* changes to an existing function (`*.m`-file), change the date and name of developer near the bottom of this commented section in the *last modified* line.

### Bugfixes, new features and functions
* For any development, whether bugfixes, introducing new functions or new/updated features for existing functions: make a separate branch from `devel` and name the branch for instance after the function/feature you are fixing/developing. If you work on a fix, start the branch name with `fix/`, if you work on a feature, start the branch name with `feat/`. Examples: `fix/format_reactions` or `feat/new_algorithms`.
* Make commits to this branch while developing. Aim for backwards compatibility to accommodate users with older software versions.
* When you are happy with your new function/feature, make a pull request to the `devel` branch. Also, see [Pull request](#pull-request) below.

### Semantic commits
Use semantic commit messages to make it easier to show what you are aiming to do:
* `chore`: updating binaries (model `MATLAB` structures), UniProt databases, physiology and protemics data files, etc.
* `doc`: updating documentation (in `doc` folder) or explanatory comments in functions.
* `feat`: new feature added, e.g. new function introduced / new parameters / new algorithm / etc.
* `fix`: bugfix.
* `refactor`: see [code refactoring](https://en.wikipedia.org/wiki/Code_refactoring).
* `style`: minor format changes of functions (spaces, semi-colons, etc., no code change).

Examples:
```
feat: exportModel additional export to YAML
chore: any routine task
fix: optimizeProb parsing results from Gurobi
```
More detailed explanation or comments can be left in the commit description.

### Pull request
* No changes should be directly commited to the `master` or `devel` branches. Commits are made to side-branches, after which pull requests are made for merging with `master` or `devel`.
* The person making the pull request and the one accepting the merge _cannot_ be the same person.
* A merge with the main branch invokes a new release.
