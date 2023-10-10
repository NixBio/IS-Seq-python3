
set -u 
set -e

# Install Homebrew firstly if you don't have brew installed on your Macs 

/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Use brew to install git and git-lfs 
brew install git
brew install git-lfs

# Under your $HOME, pull IS-Seq repository

git clone https://github.com/aiminy/IS-Seq-python3.git
cd IS-Seq-python3
git lfs pull

# Install Docker by useing the following link on your system

https://docs.docker.com/engine/install/

# Pull IS-Seq Docker image
After docker is installed, you can pull IS-Seq Docker image like the following

docker pull aiminy/isseq:1.0
