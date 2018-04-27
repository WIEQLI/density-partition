/// \file test_partition.cpp
/// \brief Test the operation of the density-based partition of sounding bin counts
///
/// This file is a simple driver to test the functionality and methods used to split up the
/// computational task of deriving a CHRT variable-resolution bathymetric estimate.

/*
 * Copyright (c) 2018, University of New Hampshire, Center for Coastal and
 * Ocean Mapping & NOAA-UNH Joint Hydrographic Center.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * (You might also be able to get a copy of the license electronically if
 * required, from <http://www.gnu.org/licenses/>)
 *
 */

#include <cstdio>
#include <string>
#include <fstream>
#include <cmath>
#include <iostream>

#include "boost/program_options.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"

#include "stdtypes.h"
#include "DensityPartition.h"

namespace po = boost::program_options;

#define VERSION_STRING	"1.0"

/// \struct Configuration
/// \brief Simple object to hold the results of a partition for a single processor count

struct Configuration {
	std::list<IntBBox>* partition;			///< Bounding boxes for the leaves of the best partition of this size
    f64                 mean_compute_time;  ///< Time to compute the partition (s)
    f64                 stddev_compute_time;///< Standard deviation of the time to compute the partition (s)
	u32					partition_cost;		///< Cost of the 'best' partition
};

void Syntax(po::options_description const& cmdopt)
{
	std::cout << "test_partition [" << VERSION_STRING << ", " << __DATE__ << ", " << __TIME__ << "] - Test density-based partitioning."
		<< std::endl;
	std::cout << "Syntax: test_partition [opt] <count_file><configs><stats><min_proc><max_proc>" << std::endl;
	std::cout << " Input count of soundings/bin -----^         ^       ^       ^         ^" << std::endl;
	std::cout << " Output file for best configurations --------'       |       |         |" << std::endl;
	std::cout << " Output file for configuration statistics -----------'       |         |" << std::endl;
	std::cout << " Minimum number of processors to consider -------------------'         |" << std::endl;
	std::cout << " Maximum number of processors to consider -----------------------------'" << std::endl;
	std::cout << cmdopt << std::endl;
}

void DumpConfigStats(std::string const& filename, std::vector<Configuration>& configs, u32 const total_observations)
{
	std::ofstream	f(filename.c_str());
	
	f << "# Columns are:" << std::endl;
	f << "# Mean Compute, Std. Dev. Compute, Partition Cost, Total Observations, Mean Observations/Thread, Std. Dev. Observations/Thread,"
		<< " Biggest Thread Observation Count, Effective Speedup, Nominal Efficiency" << std::endl;

	for (std::vector<Configuration>::const_iterator conf = configs.begin(); conf != configs.end(); ++conf) {
		u32 n_proc = static_cast<u32>(conf->partition->size());
		u32 max_thread_count = 0;	// Maximum number of soundings assigned to any one thread
		f64 mean_count, stddev_count, max_speedup;
		u32 total_idle;
		
		for (std::list<IntBBox>::const_iterator p = conf->partition->begin(); p != conf->partition->end(); ++p) {
			if (p->observation_count > max_thread_count) max_thread_count = p->observation_count;
		}
		if (n_proc > 1) {
			mean_count = static_cast<f64>(total_observations)/n_proc;
			max_speedup = static_cast<f64>(total_observations)/max_thread_count;
			total_idle = 0;
			f64 diff, sum_sq_diff = 0.0;
			for (std::list<IntBBox>::const_iterator p = conf->partition->begin(); p != conf->partition->end(); ++p) {
				diff = p->observation_count - mean_count;
				sum_sq_diff += diff*diff;
				if (p->observation_count < max_thread_count) total_idle += max_thread_count - p->observation_count;
			}
			stddev_count = sqrt(sum_sq_diff/(n_proc-1));
		} else {
			mean_count = total_observations;
			stddev_count = 0.0;
			max_speedup = 1.0;
			total_idle = 0;
		}
		f64 effective_speedup = static_cast<f64>(total_observations) / max_thread_count;
		f64 efficiency = static_cast<f64>(total_observations) / (n_proc * max_thread_count);
		f	<< conf->mean_compute_time << " " << conf->stddev_compute_time << " " << conf->partition_cost << " " << total_observations << " "
            << mean_count << " " << stddev_count << " " << max_thread_count << " " << effective_speedup << " " << efficiency << std::endl;
	}
}

void DumpConfigs(std::string const& filename, std::vector<Configuration>& configs, u32 const total_observations)
{
	std::ofstream	f(filename.c_str());
	
	for (std::vector<Configuration>::const_iterator conf = configs.begin(); conf != configs.end(); ++conf) {
		u32 n_proc = static_cast<u32>(conf->partition->size());
		u32 max_section_observations = 0;
		
		f << "--------------" << std::endl;
		f << "Number of Processors: " << n_proc << std::endl;
		for (std::list<IntBBox>::const_iterator p = conf->partition->begin(); p != conf->partition->end(); ++p) {
			f << "\t" << *p << std::endl;
			if (p->observation_count > max_section_observations)
				max_section_observations = p->observation_count;
		}
		f64 efficiency = static_cast<f64>(total_observations) / (n_proc * max_section_observations);
		f << "Nominal Efficiency: " << efficiency << std::endl;
	}
}

void DumpMachineReadable(std::string const& filename, std::vector<Configuration>& configs)
{
	std::ofstream f(filename.c_str());
	
	for (std::vector<Configuration>::const_iterator conf = configs.begin(); conf != configs.end(); ++conf) {
		u32 n_proc = static_cast<u32>( conf->partition->size() );
		
		f << n_proc << std::endl;
		for (std::list<IntBBox>::const_iterator p = conf->partition->begin(); p != conf->partition->end(); ++p) {
			f << p->left << " " << p->right << " " << p->bottom << " " << p->top << std::endl;
		}
	}
}

void DumpAllPartitions(std::ofstream& f, std::vector<Partition>& partitions, u32 const total_observations)
{	
	f << "------" << std::endl;
	f << "Total Partitions = " << partitions.size() << std::endl;
	f << "Number of Sections = " << partitions.front().size() << std::endl;
	
	u32 n_sections = partitions.front().size();
	f64 nominal_observations = static_cast<f64>(total_observations) / n_sections;
	
	for (std::vector<Partition>::const_iterator part = partitions.begin(); part != partitions.end(); ++part) {
		f << "Aspect Cost = " << part->AspectCost() << std::endl;
		f << "Count Cost = " << part->CountCost() << std::endl;
		for (std::list<DensityTree const*>::const_iterator leaf = part->begin(); leaf != part->end(); ++leaf) {
            f64 partial_cost = std::abs<f64>((*leaf)->BoundingBox().observation_count - total_observations);
			f << "\t" << "[" << *leaf << "] " << (*leaf)->BoundingBox() << "[Increment Count Cost " << partial_cost << "]" << std::endl;
		}
	}
}

int main(int argc, char **argv)
{
	po::options_description cmdopt("Options");
	cmdopt.add_options()
		("help,h",			"Generate syntax list")
		("partitions,p",	po::value<std::string>(),	"Dump file for exhaustive partition configurations")
		("tree,t",			po::value<std::string>(),	"Dump the partition tree for debug")
		("config,c",		po::value<std::string>(),	"Generate a machine-readable form of the configuration for analysis")
		("input,i",			po::value<std::string>(),	"Input file of sounding counts/bin")
		("output,o",		po::value<std::string>(),	"Output file of best partition configurations")
		("stats,s",			po::value<std::string>(),	"Output file for configuration statistics")
		("min-part,b",		po::value<int>(),			"Minimum number of partitions to consider")
		("max-part,e",		po::value<int>(),			"Maximum number of partitions to consider")
        ("binary",                                      "Read binary counts from input, rather than text")
        ("repeat,r",        po::value<int>(),           "Repeat computation and estimate compute time statistics")
		;
	po::positional_options_description cmdline;
	cmdline.add("input", 1).add("output", 1).add("stats", 1).add("min-part", 1).add("max-part", 1);
	
	po::variables_map optvals;
	po::store(po::command_line_parser(argc, argv).options(cmdopt).positional(cmdline).run(), optvals);
	po::notify(optvals);
	
	// Look for help request
	if (optvals.count("help")) {
		Syntax(cmdopt);
		return 1;
	}
	
	// Check for mandatory parameters, and convert
	if (!optvals.count("input") || !optvals.count("output") || !optvals.count("max-part")) {
        std::cout << "error: not enough mandatory parameters." << std::endl;
		Syntax(cmdopt);
		return 1;
	}
	
	// Configure the maximum number of processors to consider
	s32 min_processors = optvals["min-part"].as<s32>();
	s32 max_processors = optvals["max-part"].as<s32>();
	if (min_processors < 1 || min_processors > 100) {
		std::cout << "error: minimum number of partitions must be in the range [1,100]." << std::endl;
		return 1;
	}
	if (max_processors < 1 || max_processors > 100) {
		std::cout << "error: maximum number of partitions must be in range [1,100]." << std::endl;
		return 1;
	}
	if (max_processors < min_processors) {
		std::cout << "error: maximum number of partitions must be as least the minimum number!" << std::endl;
		return 1;
	}
	
	// Get set up for the partition dump table file if required
	std::ofstream *partition_table = NULL;
	if (optvals.count("partitions")) {
		partition_table = new std::ofstream(optvals["partitions"].as<std::string>().c_str());
	}
	// Get set up for the tree dump table file if required
	std::ofstream *tree_table = NULL;
	if (optvals.count("tree")) {
		tree_table = new std::ofstream(optvals["tree"].as<std::string>().c_str());
	}
	
	// Load input file (by the simplest means possible)
    u32 *counts;
    u32 total_observations = 0;
    u32 rows, cols;
    std::ifstream f(optvals["input"].as<std::string>().c_str());
    if (optvals.count("binary")) {
        f.read((char*)&rows, sizeof(u32));
        f.read((char*)&cols, sizeof(u32));
        std::cout << "info: Input count grid with " << rows << " rows and " << cols << " cols." << std::endl;
        counts = new u32[rows*cols];
        for (u32 n = 0; n < rows*cols; ++n) {
            f.read((char*)&counts[n], sizeof(u32));
            total_observations += counts[n];
        }
    } else {
        f >> rows >> cols;
        std::cout << "info: Input count grid with " << rows << " rows and " << cols << " cols." << std::endl;
        
        counts = new u32[rows*cols];
        for (u32 n = 0; n < rows*cols; ++n) {
            f >> counts[n];
            total_observations += counts[n];
        }
    }
	std::cout << "info: Total " << total_observations << " observations represented in counts grid." << std::endl;
	
	// Compute the configurations for processors in range [1, max_processors]
	DensityPartition partition(counts, cols, rows);
	std::vector<Configuration>	config_set;
    
    u32 n_repeats = 1;
    if (optvals.count("repeat")) {
        n_repeats = optvals["repeat"].as<int>();
        if (n_repeats < 1) {
            std::cout << "error: Number of repeats must be >= 1" << std::endl;
            return 1;
        }
    }
	
	for (s32 n_proc = min_processors; n_proc <= max_processors; ++n_proc) {
		std::cout << "info: Computing configuration for " << n_proc << " processors with " << n_repeats << " repeats." << std::endl;
		
		Configuration config;
		
        f64 sum_compute_times = 0.0, sum_sq_compute_times = 0.0;
        for (u32 r = 0; r < n_repeats; ++r) {
            boost::posix_time::ptime start(boost::posix_time::microsec_clock::universal_time());
            partition.Repartition(n_proc);
            config.partition = partition.BestPartition(config.partition_cost);
            boost::posix_time::ptime end(boost::posix_time::microsec_clock::universal_time());
            
            boost::posix_time::time_duration delay = end - start;
            sum_compute_times += delay.total_microseconds()/1.0e6;
            sum_sq_compute_times += (delay.total_microseconds()/1.0e6)*(delay.total_microseconds()/1.0e6);
            
            if (r != n_repeats-1)
                delete config.partition;
        }
        
        config.mean_compute_time = sum_compute_times/n_repeats;
        if (n_repeats > 1) {
            config.stddev_compute_time = sum_sq_compute_times/(n_repeats-1) -
                                    (static_cast<f64>(n_repeats)/(n_repeats-1))*config.mean_compute_time*config.mean_compute_time;
        } else {
            config.stddev_compute_time = -1.0;
        }
        
		config_set.push_back(config);
		
		if (partition_table != nullptr) {
			std::vector<Partition> part_list;
			partition.EnumeratePartitions(part_list);
			DumpAllPartitions(*partition_table, part_list, total_observations);
		}
		if (tree_table != nullptr) {
			partition.DumpTree(*tree_table);
		}
	}
	std::cout << "info: Complete, dumping configuration statistics." << std::endl;
	DumpConfigs(optvals["output"].as<std::string>(), config_set, total_observations);
	DumpConfigStats(optvals["stats"].as<std::string>(), config_set, total_observations);
	if (optvals.count("config"))
		DumpMachineReadable(optvals["config"].as<std::string>(), config_set);
	if (partition_table != nullptr) delete partition_table;
	if (tree_table != nullptr) delete tree_table;
	return 0;
}
