<p align="center">
    <img 
        src="https://github.com/davidmeijer/versalign/blob/main/logo.png" 
        height="150"
    />
</p>

<p align="center">
    <a href="https://github.com/davidmeijer/versalign/actions/workflows/tests.yml">
        <img 
            alt="Tests" 
            src="https://github.com/davidmeijer/versalign/actions/workflows/tests.yml/badge.svg" 
        />
    </a>
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
</p>

Versalign is Python package that allows you to create multiple sequence alignments for arbitrary lists of objects.

## 💪 Getting Started

Pairwise alignment:

```python
from versalign.motif import Motif
from versalign.pairwise import PairwiseAlignment, align_pairwise
from versalign.sequence import Sequence

class A(Motif):
    def __eq__(self, other):
        return isinstance(other, A)
    
    def __str__(self):
        return "A"

class B(Motif):
    def __eq__(self, other):
        return isinstance(other, B)
    
    def __str__(self):
        return "B"

def score_func(a, b):
    if a == b:
        return 1
    return -1

seq_a = Sequence("seq_a", [A(), A(), A()])
seq_b = Sequence("seq_b", [B(), B(), B()])

aligned_seq_a, aligned_seq_b, score = align_pairwise(
    seq_a=seq_a,
    seq_b=seq_b,
    gap_penalty=2,
    end_gap_penalty=1,
    score_func=score_func,
    algorithm=PairwiseAlignment.NEEDLEMAN_WUNSCH
)

print(aligned_seq_a)
print(aligned_seq_b)

>> AAA---
>> ---BBB
```

Multiple sequence alignment:

```python
from versalign.msa import multiple_sequence_alignment

seq_a = Sequence("seq_a", [A(), A(), A()])
seq_b = Sequence("seq_b", [B(), B(), B()])
seq_c = Sequence("seq_c", [A(), B(), B()])

result = multiple_sequence_alignment(
    seqs=[seq_a, seq_b, seq_c],
    gap_penalty=2,
    end_gap_penalty=1,
    score_func=score_func,
)

for seq in result:
    print(seq)

>> ---BBB
>> --ABB-
>> AAA---
```

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