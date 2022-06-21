# TransdEM

**TransdEM** is a package for transdimensional Bayesian inversion of electromagnetic data in layered media written in the [Julia language](http://julialang.org).

For the details regarding the algorithm and implementation, please refer to: 
> Peng, R., B. Han, Y. Liuand X. Hu, 2022, A Julia software package for transdimensional Bayesian inversion of
electromagnetic data over horizontally stratified media: Geophysics, 87(5), F29-F40; DOI:
> [10.1190/GEO2021-0534.1](https://library.seg.org/doi/10.1190/geo2021-0534.1).

*  Authors: [Ronghua Peng](https://github.com/prhjiajie), [Bo Han](https://github.com/hanbo1735), Yajun Liu and Xiangyun Hu (China University of Geosciences (Wuhan)).


## License

The **TransdEM** package is licensed under the [GNU General Public License](http://www.gnu.org/licenses/).


## File structure
* **./doc** :        an instruction of file format.

* **./examples** :   contains subdirectories corresponding to different types of synthetic examples, including all those presented in the manuscript.

* **./src** :        source code.

## Installation of Julia
TransdEM is compatible with Julia v0.7 and later versions. We recommend to install v1.0.5, the long-term support (LTS) release.

### Windows systems
Go to the [Julia download page](https://julialang.org/downloads/) to download the Windows command line version (.exe) and install it.

### Linux systems
Although Julia is a cross-platform language, we strongly recommend to run **TransdEM** under Linux rather than Windows. This is because some of the third-party packages utilized by TransdEM such as **Dipole1D** are more straightforward to complie under Linux.

There are three ways to install Julia on Linux:

* **Using precompiled binaries (recommended)**.	Go to the [Julia download page](https://julialang.org/downloads/) to download the generic Linux binaries (.tar.gz file).
	Then make sure that the Julia executable is visible for your system. To do this, first extract the .tar.gz file to a folder on your computer.
	Then you can either add Julia’s bin folder to your system PATH environment variable, or create a symbolic link to julia inside a folder which
	is on your system PATH, for example, by using the following command:

  `sudo ln -s <where you extracted the julia archive>/bin/julia /usr/local/bin/julia`


* **Compiling from source**. Assume that Git has been installed already, then we can grab the Julia sources from GitHub by using the following command:

  `git clone git://github.com/JuliaLang/julia.git`

  This will download the Julia source code into a julia directory in the current folder. The Julia building process needs the GNU compilation tools g++, gfortran, and m4, so make sure that you have installed them. Now go to the Julia folder and start the compilation process as follows:

  `cd julia`

  `make`


* **Using PPA for Ubuntu Linux**. Particularly, for Ubuntu systems (Version 12.04 or later), there is a Personal Package Archive (PPA) for Julia
	that makes the installation painless. All you need to do to get the stable version is to issue the following commands in a terminal session:

  `sudo add-apt-repository ppa:staticfloat/juliareleases`

  `sudo add-apt-repository ppa:staticfloat/julia-deps`

  `sudo apt-get update`

  `sudo apt-get install julia`

After a successful installation, Julia can be started by double-clicking the Julia executable (on Windows) or typing `julia` from the command line (on Linux). Following is the Julia's command line environment (the so-called REPL):


```jl
   _       _ _(_)_     |
  (_)     | (_) (_)    |  Documentation: https://docs.julialang.org
  _  _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.0.5 (2019-09-09)
 _/ |\__'_|_|_|\__'_|  |  Official http://julialang.org/ release
|__/                   |  

julia>
```

## Installation of the TransdEM package
### Setting up the package environment
**TransdEM** depends on several external packages (the so-called dependencies) which are not shipped with the package .zip file. These dependencies can be automatically resolved by activating and instantiating the package environment through [Julia's package manager (Pkg)](https://julialang.github.io/Pkg.jl/v1/). Go to the package directory, for example:  
`cd /home/username/code/TransdEM`

, and then enter the Julia REPL. Then press `]` from the Julia REPL you will enter the Pkg REPL which looks like
```jl
(v1.0) pkg>
```

, indicating that you are currently in the environment named v1.0, Julia 1.0's default environment. To switch to the package environment, just `activate` the current directory:
```jl
(v1.0) pkg> activate .
```
you will get:
```jl
(TransdEM) pkg>
```
indicating that you are in the environment TransdEM. The environment will not be well-configured until you `instantiate` it:
```jl
(TransdEM) pkg> instantiate
```
. By doing so the dependencies listed in `Project.toml` and `Manifest.toml` can be automatically downloaded and installed.

To get back to Julia REPL from Pkg REPL, press `backspace` or `^C`.

### Building the em1dmod library
**TransdEM** contains several Fortran codes that are used to compute 1D MT/CSEM/TEM responses, and they must be compiled to generate a shared library (.so, on Linux) or dynamic link library (.dll, on Windows) called `em1dmod` prior to running **TransdEM**. Suppose that you have already installed the GNU fortran compiler (gfortran) in your computer system, all you need to do is going to the subdirectory of the `TBFwdSolver` module and "including" the script *build.jl*, the Julia REPL commands of which are like:
  ```jl
  julia> cd("/home/username/code/TransdEM/src/TBFwdSolver/deps")
  julia> include("build.jl")
  ```
  If no error is reported during the building process, the library file should be generated and placed at the directory like “D:\code\TransdEM\src\TBFwdSolver\deps\lib”. 
  
  For installing gfortran, on Windows systems we recommend to install [MinGW](https://en.wikipedia.org/wiki/MinGW/), particularly [Mingw-w64](https://www.mingw-w64.org/) for 64 bits systems. The windows installer (.exe) can be downloaded from https://sourceforge.net/projects/mingw-w64. After installing, you need to edit the `PATH` variable. You can access the System Control Center by pressing **Windows Key + Pause**. In the System window, click **Advanced System Settings $\rightarrow$ Advanced (tab) $\rightarrow$ Environment Variables**. For Windows 10, a quick access is to enter "Edit the system environment variables" in the Start Search of Windows and click the button "Environment Variables". Change the `PATH` variable (double-click on it or Select and **Edit**), and add the path where your MinGW-w64 has been installed to e.g., `C:\mingw\mingw64\bin`. This folder should contain a number of .exe-files that you can see in your explorer. To check that your Mingw-w64-gfortran is correctly installed and available, enter the shell mode by pressing ";" from the Julia REPL and type `gfortran --version`, the expected output looks like
```jl
shell> gfortran --version
GNU Fortran (x86_64-posix-seh-rev0, Built by MinGW-W64 project) 8.1.0
Copyright (C) 2018 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
```
## Running the code
### Running single MCMC sampling
* First, you need to let the **TransdEM** package to be "loaded" by the current Julia environment. This is done by adding the parent directory of the package directory to  `LOAD_PATH`, a global environment variable of Julia. For example, the TransdEM package is placed at `home/username/code` on Linux or at `D:\\code` on Windows, then type the following command from the Julia REPL:

  ```jl
  julia> push!(LOAD_PATH,"/home/username/code")
  ```

  on Linux, or

  ```jl
  julia> push!(LOAD_PATH,"D:\\code")
  ```

  on Windows.   


* Then go to the directory where the running script loaded, and run the script by typing the following command (for example) from the Julia REPL:

  ```jl
  julia> include("runMCMCScript.jl")
  ```

### Running parallel MCMC sampling
To perform parallel MCMC sampling, call the parallel MCMC sampling function **parallelMCMCsampling** instead of  the single MCMC sampling **runMCMC** (please refer to the `paraMCMCScript.jl` scripts within the `examples` directory).

* first you need to launch multiple worker processes by either starting Julia like

  `shell> julia -p 12`

  or adding processes within Julia (recommended) like

  `julia> addprocs(12)`

* Then you need to get the **TransdEM** package to be "loaded" on all processes by following command from the Julia REPL:

	```jl
	julia> @everywhere push!(LOAD_PATH,"/home/username/code")
	````

* Finally, go to the directory where the running script loaded, and run the script by typing the following command (for example) from the Julia REPL:

  ```jl
  julia> include("paraMCMCScript.jl")
  ```

### Running MCMC sampling with parallel tempering
In order to accelerate convergence of the Markov chains, a powerful technique known as parallel tempering can be used. To perform MCMC sampling with parallel tempering, call the parallel tempered MCMC sampling function **runTemperedMCMC** instead of the parallel MCMC sampling function **parallelMCMCsampling** (please refer to the `runPTMCMCScript.jl` scripts within the `examples` directory).
* first you need to launch multiple worker processes by either starting Julia like

  `shell> julia -p 6`

  or adding processes within Julia (recommended) like

  `julia> addprocs(6)`
Note that the number of processes invoked should be equal or larger than the number of parallel tempered Markov chains.
* Then you need to get the **TransdEM** package to be "loaded" on all processes by following command from the Julia REPL:

	```jl
	julia> @everywhere push!(LOAD_PATH,"/home/username/code")
	````

* Finally, go to the directory where the running script loaded, and run the script by typing the following command (for example) from the Julia REPL:

  ```jl
  julia> include("runPTMCMCScript.jl")
  ```

### Writing a running script
For each numerical example contained in the directory `./examples`, one or more running scripts named `runMCMCScript.jl`, `paraMCMCScript.jl` or `runPTMCMCScript.jl` have been provided. These scripts are well documented. A user can modify them to get his/her own.
