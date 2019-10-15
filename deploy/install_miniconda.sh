# Install miniconda
MINICONDA_URL="http://repo.anaconda.com/miniconda"

if [[ $TRAVIS_OS_NAME == "linux" ]]; then
    MINICONDA_FILE="Miniconda3-latest-Linux-x86_64.sh"
else
    MINICONDA_FILE="Miniconda3-latest-MacOSX-x86_64.sh"
fi

wget $MINICONDA_URL/$MINICONDA_FILE -O $HOME/miniconda.sh
bash $HOME/miniconda.sh -b -p $HOME/miniconda

# Add miniconda to PATH
export PATH=$HOME/miniconda:$PATH
export PATH=$HOME/miniconda/Scripts:$PATH
export PATH=$HOME/miniconda/bin:$PATH
