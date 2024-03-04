# maligner: An Open-Source Molecular Alignment Tool

![GitHub top language](https://img.shields.io/github/languages/top/hellmrf/maligner)
[![GitHub License](https://img.shields.io/github/license/hellmrf/maligner)](https://github.com/hellmrf/maligner/blob/main/LICENSE)
[![Powered by RDKit](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)


## Current plan

Currently, the plan for the UX of the software is the following.

**Main screen**:
- The user can load molecules (in .mol, .mol2 or ??)
- Can see a grid with all molecules and the filename
- Can mark one molecule as anchor
- Can go to structure selector.
- Options menu
- "Preview Alignment" and "Save Alignment" buttons.

**Options menu**
- Save action:
  - Overwrite original files
  - Suffix aligned files (option to write the suffix)
  - *(Default)* Save in a separate directory (default: `aligned`)

**Atom selector**:
- User see the molecule and, by default, the MCS calculated.
- Three selection buttons: Select MCS, Select All, Select None.
- Two more buttons: Cancel, Save.
- User can click on atoms to select/deselect.
- Select structure by SMILES/SMARTS.
- Maybe implement lasso tool in future.

