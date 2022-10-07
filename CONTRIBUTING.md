# Contributing

## Setting up a development environment

- Install development dependencies.

```
python3 -m pip install -r requirements-dev.txt
```

- Install [pre-commit](https://pre-commit.com/) hooks.

```
python3 -m pre_commit install
```

## Running Pylint

Use [Pylint](https://www.pylint.org/) to check for some types of errors.

```
./lint
```

To disable warnings, use:

```
./lint --disable=W
```

## Running Black

Use [Black](https://black.readthedocs.io/) to format code.

```
black gnomad
```

## Running isort

Use [isort](https://pycqa.github.io/isort/) to format imports.

```
isort gnomad
```

## Running pydocstyle

Use [pydocstyle](https://www.pydocstyle.org/) to check that docstrings conform to [PEP 257](https://www.python.org/dev/peps/pep-0257/).

```
pydocstyle gnomad
```
