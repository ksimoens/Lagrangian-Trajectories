# Setup netCDF Library for C++

Instructions on how to set up the netCDF library for Ubuntu (version 18.04 and above).  

## Ubuntu netCDF Library

If not already installed, first install the general netCDF library:  
```
$ sudo apt-get libnetcdf-dev
```  
Check whether the netCDF library has been correctly installed by:  
```
$ nc-config --has-nc4
```  
which should return:  
```
yes
```

## NetCDF Library for C++

1. Download the desired library [release](https://downloads.unidata.ucar.edu/netcdf/) and unzip the directory.  
Setup instructions can be found on the [GitHub](https://github.com/Unidata/netcdf-cxx4/) page.

2. Here, we will follow the `cmake` instructions.  
If not yet installed, install `cmake` by:  
```
$ sudo apt-get install cmake
``` 

3. Set up the Cmake Makefile. In the `netcdf-cxx4-X.X.X` directory, create:  
```
$ mkdir build  
$ cd build
```  
The Cmake Makefile is created by the Cmake command:  
```
$ cmake ..
```  
However, it is possible that the required `hdf5` library is not found.  
First of all, make sure that the library is properly installed:  
```
$ sudo apt-get install libhdf5-serial-dev
```  
To check the hdf5 version, use:  
```
$ dpkg -l | grep hdf5
```  
However, it is possible that the library is still not found by Cmake because the INCLUDE variables are not automatically set correctly. To remedy this, run the `cmake` command setting the INCLUDE variables manually:  
```
$ cmake .. -DCMAKE_C_FLAGS="-I /usr/include/hdf5/serial/" -DCMAKE_PREFIX_PATH=/usr
```  
where the default library locations are used. Change the locations to the hdf5 library locations on your distribution.

4. Make the netCDF C++ library. From within the same `build` directory, build the library:  
```
$ make  
$ ctest  
$ make install
```  
Finally, in order to be able to link the netCDF library in your C++ script, run the [ldconfig](https://www.quora.com/What-does-ldconfig-do) command:  
```
$ sudo ldconfig
```

## Compilation with netCDF C++ Library

In order to use the netCDF C++ interface, the `lnetcdf-cxx4` and `lnetcdf` libraries have to be linked at compile time. Additionally, the location of these libraries has to be specified, for example `/usr/local/lib` is the  default location of these libraries, but could be different on your distribution.  
A compilation example could then be:  
```
$ g++ main.cpp my_netcdf_script.cpp -L/usr/local/lib -lnetcdf-cxx4 -lnetcdf -o main
```  
Further examples can be found in the library directory or in the Makefile script of this programme.