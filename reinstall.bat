pip uninstall maligner -y
hatch build -t sdist
pip install dist/maligner-0.0.0.tar.gz