[![CI](https://github.com/schrodinger/pymol-open-source/workflows/CI/badge.svg)](https://github.com/schrodinger/pymol-open-source/actions)

<img src="./data/pymol/icons/icon2.svg" height="100" align="right" />

# Open-Source PyMOL

[Open-source foundation](https://pymol.org/#opensource) of the user-sponsored PyMOL molecular visualization system.

The commercial PyMOL product ("Incentive PyMOL") with maintenance and support is available from https://pymol.org

## Installation

See [INSTALL](INSTALL).

#### TL,DR;

Clone repo and change to the directory
```
git clone https://github.com/schrodinger/pymol-open-source
cd pymol-open-source/
```

Install with base dependencies
```
sudo apt update
sudo apt install cmake netcdf-bin libnetcdf-dev libglm-dev libglew-dev libpng-dev libfreetype-dev libfreetype6 libfreetype6-dev
pip install .
```

Install Qt5
```
sudo apt install qtcreator qtbase5-dev qt5-qmake cmake pyqt5-dev pyqt5-dev-tools python3-pyqt5
pip install pyqt5
pymol
```
![image](https://github.com/user-attachments/assets/b9388abb-d525-4ced-87fd-6396d1f44152)

## Contributing

See [DEVELOPERS](DEVELOPERS).

## License

Copyright (c) [Schrodinger, LLC](https://www.schrodinger.com/)

Published under a BSD-like license, see [LICENSE](LICENSE).
