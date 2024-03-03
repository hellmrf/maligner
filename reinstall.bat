pip uninstall malign -y
hatch build -t sdist
pip install dist/malign-0.0.0.tar.gz