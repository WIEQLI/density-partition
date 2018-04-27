Center for Coastal and Ocean Mapping & NOAA/UNH Joint Hydrographic Center
University of New Hampshire
Durham, NH 03820

April 25, 2018

This example code distribution accompanies the paper "Parallel Implementation of a Variable Resolution Bathymetric Estimator", B.R. Calder, published in Computers and Geosciences.  The code demonstrates the ideas expressed in the paper on partitioning the computational problem attempted into a (small) number of segments that have approximately the same amount of work.  This allows the code to throw off a thread for each segment without having to coordinate the edges with the others, or continually communicate with the controller thread for more activity in order to load-balance the computation.  Note that the actual computational code is not provided (see "Computationally Efficient Variable Resolution Depth Estimation", B.R. Calder and G. Rice, Computers and Geosciences 106:49-59, 2017 for more details, and example code).

The code is organised in a primary C++ class (DensityPartition) and a collection of helper classes to implement the partitioning of the whole computation problem with respect to data density.  The test_partition.cpp code provides a simple command-line driver to read a data density estimate grid, and run the partitioning process with control over the number of partitions to generate, and options for generating debugging information and making machine readable versions of the output so that the statistics of the partition can be further analysed.  Timing code (with statistics gathering) is also provided.

This code is copyright of the University of New Hampshire as indicated in the source files, and is released under version two of the GNU General Public License.  Please see LICENSE.txt for the full text of the license.

File Descriptions
-------------------

The code consists of a set of classes in DensityPartition.cpp which implement the partitioning algorithm, with corresponding declarations in DensityPartition.h.  Standard types for simple objects with known size are provided in stdtypes.h.  All of the code is extensively documented using Doxygen conventions.

The code in test_partition.cpp is used to drive code from the command line in order to examine the algorithm and provide debugging information to check on the behaviour of the code.  In addition to exercising the algorithm over a number of different partition counts (typically one for each thread of execution to be attempted), the code can repeat the computation for each partition count, and accumulate statistics on the time taken.

All of the various files generated are simple ASCII, and are mostly self-explanatory.  The output from the "--config" option, however, is intentionally kept simple (i.e., without human-readable documentation) so that it is easier to parse in code for further analysis/plotting/etc.  The format is simple, and can be found in the DumpMachineReadable() free-function in test_partition.cpp.  This is, however, equivalent to the information provided through the DumpConfigs() free-function, which is the default for the "--output" option.

Compiling the Code
-----------------------

The code is simple C++11, and should be able to be compiled on most modern C++ compilers.  The primary development system was macOS using LLVM 8.0.0 (clang-800.0.42.1), although the code is known to also compile and run on Windows and Linux.  The example driver code uses Boost for programme options handling and timestamp management; these dependencies could be removed fairly easily if required, however.

Building the code is simply a matter of compiling and linking.  Assuming an installation of the Boost library in ${HOME}/include/boost and ${HOME}/lib, the command-line used to compile the code will be typically:

c++ -o test_partition -I~/include test_partition.cpp DensityPartition.cpp -L~/lib -lboost_date_time -lboost_program_options

Test Data and Results
--------------------------

A pair of example datasets, corresponding to the surveys described in the paper, are provided for testing the code, along with the expected outputs.  Test data is provided in ./TestData, and the corresponding results in ./TestResults.  The input data files are binary encoded counts of observations in the low-resolution grid cells computed in the first phase of the CHRT algorithm to determine the data density.  The command-line options for the test runs were:

../test_partition 	--partitions h11825_all_partitions.txt
				--tree h11825_partition_tree.txt
				--config h11825_configuration.txt
				--input ../TestData/h11825_counts.raw
				--output h11825_best_partition.txt
				--stats h11825_stats.txt
				--min-part 1
				--max-part 8
				--binary
				--repeat 100

../test_partition	--partitions h11077_all_partitions.txt
				--tree h11077_partition_tree.txt
				--config h11077_configuration.txt
				--input ../TestData/h11077_counts.raw
				--output h11077_best_partition.txt
				--stats h11077_stats.txt
				--min-part 1
				--max-part 8
				--binary
				--repeat 100