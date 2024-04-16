# Smilei Installation Guide

The supported platforms are Linux and MacOS. If you use Windows, follow these instructions to set up a running Linux platform:

- Go to Start. Search for "Turn Windows features on or off."
- Check the option Windows Subsystem for Linux.

  <img src="https://github.com/CompuNanoLab/PlasmaPhysicsLab/assets/140382467/75fe59a8-35c6-47c9-9b3e-fef7ff2c0fad" width="400" />
- Open Command Prompt as an administrator.

  <img src="https://github.com/CompuNanoLab/PlasmaPhysicsLab/assets/140382467/96052034-10c0-417a-90dc-71180df4704d" width="400" />
  
- Run the command below to install the Distro (e.g. Ubuntu, Debian, ...) of your choice:
  ```bash
  wsl --install -d <Distro>
  ```
- Launch the Distro by searching from the start menu and insert a username and password.

## General dependencies
First, install on your laptop using a package manager the following fundamental software (dependencies) necessary to compile and run Smilei:

### Linux

For **Debian-based (Ubuntu) OS**:
```bash
sudo apt-get update
```
```bash
sudo apt-get install git python3-h5py python3-ipython python3-pint python3-sphinx python3-matplotlib python3-dev python3-numpy python3-scipy python3-pip build-essential gcc libhdf5-openmpi-dev 
```
and add the following lines to your ~/.bashrc or ~/.bash_profile file
```
export PYTHONEXE=python3
export HDF5_ROOT_DIR=/usr/lib/x86_64-linux-gnu/hdf5/openmpi
```

For **Fedora OS**:
```bash
sudo dnf install gcc-c++ git hdf5-openmpi hdf5-openmpi-devel openmpi-devel python python-devel python3-h5py ipython python3-pint python3-sphinx python3-matplotlib 
```
and add the following lines to your ``~/.bashrc`` or ``~/.bash_profile`` file

```bash
module load mpi
```

 For **ArchLinux OS**:
```bash
sudo pacman -S git hdf5-openmpi python-numpy python-sphinx python-h5py-openmpi python-matplotlib python-pint make gcc 
```

If you encounter problems, you may need to install openmpi and/or hdf5 directly from source. To do that, try to follow the instructions [here](https://smileipic.github.io/Smilei/Use/install_linux.html#troubleshooting).
 
### MacOS

First install Homebrew via:
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
Once installed, to use Homebrew on the command line you need to modify the ".zprofile" on your home by running the following commands:
```bash
(echo; echo 'eval "$(/opt/homebrew/bin/brew shellenv)"') >> /Users/<your_account>/.zprofile
eval "$(/opt/homebrew/bin/brew shellenv)""
```
Then install git and Python using Homebrew:
```bash
brew install git python 
```
To use the installed Python as the default one you will need to modify the ".zprofile" by adding the following line:
```bash
export PATH="/opt/homebrew/opt/python@3.12/libexec/bin:$PATH"
```
the path may change, use the one shown at the end of the installation. Once Python has been installed on your Laptop you have to create a virtual environment on your home to be able to install the Python packages via pip:

``` Bash
python3 -m venv myenv
```
then, if you want to, let the terminal source it automatically every time you open a new terminal by adding to your ".zshrc" file the following line:

```
source myenv/bin/activate
```
Alternatively, use the same line directly on the command line of your terminal. Source the environment and install the packages via pip:
```
pip install h5py ipython pint sphinx matplotlib dev numpy scipy 
```
Then, follow the related instructions for each case.

To install the necessary dependencies for Smilei :
```
brew update
brew install openmpi hdf5-mpi libomp adios2 ccache cmake fftw pkg-config openblas
```
The dependencies will install also *numpy*, which you should have already installed when creating the virtual environment for Python. You'd rather uninstall it by using:
```
brew uninstall --ignore-dependencies numpy
```
To check the formulae (dependencies) installed on your Mac use the command
```
brew list
```
To check the version of your C++ compiler use the command:
```
brew info gcc
```
## Build on Linux

Download the source code to the path of your choice:
```
cd /path/of/your/choice/
git clone https://github.com/SmileiPIC/Smilei.git smilei
```

Move to the downloaded directory:
```
cd smilei
```

then compile:
```
make -j 2 config="omptasks"
``` 

build the Python module to manage the output: 
```
make happi
```

Now in a Python script you can now do:
```
import happi
```

# Build on Mac

Update the profile `~/.zprofile` (or `~/.bash_profile` in some older versions of MacOS):
```
cd
nano .zprofile
```
and add the following lines:
```
export OMPI_CXX=g++-13
export HDF5_ROOT_DIR=/opt/homebrew/opt/hdf5-mpi
export PYTHONEXE=python3
```
You may need to change the `g++-13` to some other version you have on your laptop. To check the version use the `brew info gcc` command, do not forget the number of the version: by using in fact `g++` only, you will call `clang++` compiler which is provided by Apple. Such a compiler does not work with `openmpi` that is necessary for Smilei compilation instead.

Use then `git` to copy Smilei on your home:
```
git clone https://github.com/SmileiPIC/Smilei.git smilei
```
move in the folder and use `make` to copile the source (remind to activate the python environment):
```
cd smilei
make -j 8
```
To use the post-processing tools offered by Smilei, once you have compiled the source:
```
make happi
```

# Run 
To run efficiently in parallel on your machine you need to know your architecture.
For example, you can find out the number of threads per core and cores per socket on your machine with the command: `lscpu` (Linux) or `sysctl -a | grep machdep.cpu` (MacOS).

Copy from the `smilei/bin` directory the executable (smilei) and an input file (input.py) in a directory of your choice and move there.
Set the number of threads per core depending on the machine 
for example, if Thread(s) per core = 2 (in the output of the `lscpu` command), then
```
export OMP_NUM_THREADS=2
```

and if you have Core(s) per socket = 4, you can run on the 4 cores like this 
```
mpirun -np 4 ./smilei input.py
```
