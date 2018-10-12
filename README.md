# fpdb
fetch PDB database information

# install
There are several method to use fpdb
## method 1. setup.py
```bash
# > python2.7
python setup.py install
# if without root
python setup.py install --user
```
## method 2. add softlink  fpdb/fpdb to your PYTHONPATH
```bash
mkdir ~/bin
ln -s  path_to/fpdb/fpdb ~/bin/fpdb
echo "export PYTHONPATH=~/bin:$PYTHONPATH" >> ~/.bashrc
source ~/.bashrc
```
## method 3. cp fpdb/fpdb to your work dir.
