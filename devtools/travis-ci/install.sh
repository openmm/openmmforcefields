# Temporarily change directory to $HOME to install software
pushd .
cd $HOME

if [[ "$TRAVIS_OS_NAME" == "osx" ]];   then MINICONDA=Miniconda3-latest-MacOSX-x86_64.sh; fi
if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then MINICONDA=Miniconda3-latest-Linux-x86_64.sh;  fi

MINICONDA_MD5=$(curl -s https://repo.continuum.io/miniconda/ | grep -A3 $MINICONDA | sed -n '4p' | sed -n 's/ *<td>\(.*\)<\/td> */\1/p')
wget https://repo.continuum.io/miniconda/$MINICONDA
bash $MINICONDA -b
rm -f $MINICONDA

export PATH=$HOME/miniconda3/bin:$PATH

conda update -yq conda
conda install -yq conda-build jinja2 anaconda-client pip

# Restore original directory
popd
