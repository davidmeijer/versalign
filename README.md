<p align="center">
    <img 
        src="https://github.com/davidmeijer/versalign/blob/main/logo.png" 
        height="150"
    >
</p>

<p align="center">
    <a href="https://github.com/davidmeijer/versalign/actions/workflows/tests.yml">
        <img 
            alt="Tests" 
            src="https://github.com/davidmeijer/versalign/actions/workflows/tests.yml/badge.svg" 
        />
    </a>
    <!-- <a href="https://pypi.org/project/versalign">
        <img 
            alt="PyPI" 
            src="https://img.shields.io/pypi/v/versalign" 
        />
    </a> -->
    <!-- <a href="https://pypi.org/project/versalign">
        <img 
            alt="PyPI - Python Version" 
            src="https://img.shields.io/pypi/pyversions/versalign" 
        />
    </a> -->
    <!-- <a href="https://github.com/davidmeijer/versalign/blob/main/LICENSE">
        <img 
            alt="PyPI - License" 
            src="https://img.shields.io/pypi/l/versalign" 
        />
    </a> -->
    <!-- <a href="https://codecov.io/gh/davidmeijer/versalign/branch/main">
        <img 
            src="https://codecov.io/gh/davidmeijer/versalign/branch/main/graph/badge.svg" 
            alt="Codecov status" 
        />
    </a>   -->
    <a href="https://github.com/cthoyt/cookiecutter-python-package">
        <img 
            alt="Cookiecutter template from @cthoyt" 
            src="https://img.shields.io/badge/Cookiecutter-snekpack-blue" 
        />
    </a>
    <a href="https://github.com/psf/black">
        <img 
            src="https://img.shields.io/badge/Code%20style-black-000000.svg" 
            alt="Code style: black" 
        />
    </a>
    <a href="https://github.com/davidmeijer/versalign/blob/main/.github/CODE_OF_CONDUCT.md">
        <img 
            src="https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg" 
            alt="Contributor Covenant"
        />
    </a>
    <!-- <a href="https://doi.org/<doi>">
        <img 
            src="https://zenodo.org/badge/DOI/<doi>.svg" 
            alt="DOI"
        />
    </a> -->
</p>

Versalign is Python package that allows you to create multiple sequence alignments for arbitrary lists of objects.

<!-- ## 💪 Getting Started

... -->

## 🚀 Installation

<!-- The most recent release can be installed from
[PyPI](https://pypi.org/project/versalign/) with:

```shell
pip install versalign
``` -->

The most recent code and data can be installed directly from GitHub with:

```shell
pip install git+https://github.com/davidmeijer/versalign.git
```

## 👐 Contributing

Contributions, whether filing an issue, making a pull request, or forking, are appreciated. See
[CONTRIBUTING.md](https://github.com/davidmeijer/versalign/blob/main/.github/CONTRIBUTING.md) for more information on getting involved.

## 👋 Attribution

### ⚖️ License

The code in this package is licensed under the MIT License.

### 🍪 Cookiecutter

This package was created with [@audreyfeldroy](https://github.com/audreyfeldroy)'s
[cookiecutter](https://github.com/cookiecutter/cookiecutter) package using [@cthoyt](https://github.com/cthoyt)'s
[cookiecutter-snekpack](https://github.com/cthoyt/cookiecutter-snekpack) template.

## 🛠️ For Developers

<details>
  <summary>See developer instructions</summary>

The final section of the README is for if you want to get involved by making a code contribution.

### Development Installation

To install in development mode, use the following:

```bash
git clone git+https://github.com/davidmeijer/versalign.git
cd versalign
pip install -e .
```

### 🥼 Testing

After cloning the repository and installing `tox` with `pip install tox`, the unit tests in the `tests/` folder can be
run reproducibly with:

```shell
tox
```

Additionally, these tests are automatically re-run with each commit in a
[GitHub Action](https://github.com/davidmeijer/versalign/actions?query=workflow%3ATests).

### 📦 Making a Release

After installing the package in development mode and installing
`tox` with `pip install tox`, the commands for making a new release are contained within the `finish` environment
in `tox.ini`. Run the following from the shell:

```shell
tox -e finish
```

This script does the following:

1. Uses [Bump2Version](https://github.com/c4urself/bump2version) to switch the version number in the `setup.cfg`,
   `src/versalign/version.py`, and [`docs/source/conf.py`](docs/source/conf.py) to not have the `-dev` suffix
2. Packages the code in both a tar archive and a wheel using [`build`](https://github.com/pypa/build)
3. Uploads to PyPI using [`twine`](https://github.com/pypa/twine). Be sure to have a `.pypirc` file
   configured to avoid the need for manual input at this step
4. Push to GitHub. You'll need to make a release going with the commit where the version was bumped.
5. Bump the version to the next patch. If you made big changes and want to bump the version by minor, you can
   use `tox -e bumpversion -- minor` after.

</details>