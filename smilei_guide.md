# Smilei Installation Guide

The supported platforms are Linux and MacOS. If you use Windows, follow these instructions to set up a running Linux platform:

- Go to Start. Search for "Turn Windows features on or off."
- Check the option Windows Subsystem for Linux.

  <img src="https://github.com/CompuNanoLab/PlasmaPhysicsLab/assets/140382467/75fe59a8-35c6-47c9-9b3e-fef7ff2c0fad" width="400" />
- Open Command Prompt as an administrator.

  <img src="https://github.com/CompuNanoLab/PlasmaPhysicsLab/assets/140382467/96052034-10c0-417a-90dc-71180df4704d" width="400" />
  
- Run the command below to install the Distro (e.g. Ubuntu, Debian, ...) of your choice:
  ```
  wsl --install -d <Distro>
  ```
- Launch the Distro by searching from the start menu and insert a username and password.
- 
  <img src="https://github.com/CompuNanoLab/PlasmaPhysicsLab/assets/140382467/fb762f06-3d6b-4ddc-b090-e2569d73d3fa" width="400" />

## Install on Linux

### General dependencies
First, install on your laptop using a package manager the following fundamental software (dependencies) necessary to compile and run Smilei. With the same steps, you are also installing Python which will allow you to run Python scripts (.py) and jupyter notebooks.

For **Debian-based (Ubuntu) OS**:
```
sudo apt-get update
sudo apt-get install git python3-h5py python3-ipython python3-pint python3-sphinx python3-matplotlib python3-dev python3-numpy python3-scipy python3-pip build-essential gcc libhdf5-openmpi-dev 
```
Open your ``.bashrc`` or ``.bash_profile`` file in your $HOME. For example, if you want to use the nano editor, type in the terminal:
```
nano $HOME/.bashrc
```
Add the following lines:
```
export PYTHONEXE=python3
export HDF5_ROOT_DIR=/usr/lib/x86_64-linux-gnu/hdf5/openmpi
```
and source it:
```
source $HOME/.bashrc
```

For **Fedora OS**:
```
sudo dnf install gcc-c++ git hdf5-openmpi hdf5-openmpi-devel openmpi-devel python python-devel python3-h5py ipython python3-pint python3-sphinx python3-matplotlib 
```
Open your ``.bashrc`` or ``.bash_profile`` file in your $HOME. For example, if you want to use the nano editor, type in the terminal:
```
nano $HOME/.bashrc
```
Add the following lines:
```
module load mpi
export HDF5_ROOT_DIR=/usr/lib64/openmpi/
```
and source it:
```
source $HOME/.bashrc
```

 For **ArchLinux OS**:
```
sudo pacman -S git hdf5-openmpi python-numpy python-sphinx python-h5py-openmpi python-matplotlib python-pint make gcc 
```

If you encounter problems, you may need to install openmpi and/or hdf5 directly from source. To do that, try to follow the instructions [here](https://smileipic.github.io/Smilei/Use/install_linux.html#troubleshooting).

### Build

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
make -j 2
``` 

build the Python module to manage the output: 
```
make happi
```

Now in a Python script you can do:
```
import happi
```

## Install on Mac

First, you will need to install Xcode and the Command Line Tools to be able to compile Smilei:
```
xcode-select --install
```
and follow the instructions.

Here we show how to install all dependencies needed by Smilei using Brew. 

Install HomeBrew:
```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```
Install Smilei:
```
brew install --HEAD iltommi/brews/smilei
```
Smilei executables (smilei and smilei_test) and the Python module are now accessible from everywhere.

Install Python packages needed for the happi python module:
```
pip3 install ipython h5py pint sphinx matplotlib scipy
```

## Run 
To run efficiently in parallel on your machine you need to know your architecture.
For example, you can find out the number of threads per core and cores per socket on your machine with the command: `lscpu` (Linux) or `sysctl -a | grep machdep.cpu` (MacOS).

Copy from the `smilei/bin` directory the executable (smilei) and an input file (input.py) in a directory of your choice and move there.
Set the number of threads per core depending on the machine 
for example, if Thread(s) per core = 2 (in the output of the `lscpu` command), then
```
export OMP_NUM_THREADS=2
export OMP_SCHEDULE=dynamic
```

and if you have Core(s) per socket = 4, you can run on the 4 cores like this 
```
mpirun -np 4 ./smilei input.py
```
