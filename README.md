<div align="center">
  <!--<img src="docs/images/logo.svg" alt="Logo" width="100" height="100">-->

  # `maligner`
  
  _An Open-Source Molecular Alignment Tool._

  ![GitHub top language](https://img.shields.io/github/languages/top/hellmrf/maligner)
  [![GitHub License](https://img.shields.io/github/license/hellmrf/maligner)](https://github.com/hellmrf/maligner/blob/main/LICENSE)
  [![Powered by RDKit](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)
  
  <!--<a href="#about"><strong>Explore the docs »</strong></a>
  <br />
  <br />
  <a href="https://github.com/hellmrf/maligner/issues/new?assignees=&labels=bug&template=01_BUG_REPORT.md&title=bug%3A+">Report a Bug</a>
  ·
  <a href="https://github.com/hellmrf/maligner/issues/new?assignees=&labels=enhancement&template=02_FEATURE_REQUEST.md&title=feat%3A+">Request a Feature</a>-->
</div>

## About

maligner is an open-source tool for molecular alignment, designed to facilitate the manipulation and visualization of molecular structures. The goal is to provide an intuitive graphical interface and advanced functionalities to align molecules efficiently, leveraging the power of RDKit.

### Built With

- [Python](https://www.python.org/)
- [RDKit](https://www.rdkit.org/)

## Getting Started (development)

### Prerequisites

- Install [`uv`](https://docs.astral.sh/uv).

A global Python installation is not required as `uv` will handle the environment setup.

## Usage

To run the application, use the following command:

```bash
uv run gui
```

This will:
1. Create and resolve a new environment (located at `.venv`).
2. Install all necessary dependencies (including dev-dependencies).
3. Launch the application's GUI.

## Development

### Stage of development

The project is currently in the `IDEA` stage and is not yet ready for real usage. Contributions are welcome and appreciated!

### Development Setup

We use the following dev-only tools.

- [`mypy`](https://mypy-lang.org/) for type-checking.
- [`ruff`](https://github.com/charliermarsh/ruff) for linting and formatting.

Run the following commands as needed to run them.

```bash
$ uv run mypy
$ uv run ruff check
$ uv run ruff format
```

## Roadmap

- **Main Screen**:
  - [X] Load molecules in `.mol`, `.mol2`, etc.
  - [X] Display molecules in a grid with filenames.
  - [X] Mark a molecule as anchor.
  - [X] Navigate to structure selector.
  - [ ] Add options menu.
  - [ ] Implement "Preview Alignment" and "Save Alignment" buttons.

- **Options Menu**:
  - Overwrite original files.
  - Suffix aligned files (with user-defined suffix).
  - *(Default)* Save in a separate directory (`aligned`).

- **Atom Selector**:
  - [X] Display molecule and calculate MCS by default.
  - [ ] Add selection buttons: Select MCS, Select All, Select None.
  - [ ] Include Cancel and Save buttons.
  - [X] Allow users to click on atoms to select/deselect.
  - [ ] Enable selection by SMILES/SMARTS.
  - [ ] (Optional) Implement lasso tool for selection.

## Support

For support, reach out via:

- [GitHub Issues](https://github.com/hellmrf/maligner/issues)

## License

This project is licensed under the GPL-3.0 License. See the [LICENSE](https://github.com/hellmrf/maligner/blob/main/LICENSE) file for details.