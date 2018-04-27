/// \file DensityPartition.h
/// \brief Provide an object to partition sounding count arrays such that they can be computed in parallel.
///
/// This file provides an object to assist in computing the semi-optimal partition of a grid of sounding
/// counts in SuperGrid cells into N sections, such that each section contains as closely as possible the
/// same number of soundings, and consists of as close to square as possible bounding boxes.  The partition
/// splits on the SuperGrid cells, and constrains the sections to be rectangular.

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

#ifndef __DENSITY_PARTITION_H__
#define __DENSITY_PARTITION_H__

#include <list>
#include <vector>

#include "stdtypes.h"

struct CountData;

/// \struct IntBBox
/// \brief A simple integer bounding box with natural number limits
///
/// This holds the bounding box limits in SuperGrid cells for a section of the working space.  The array bounds are inclusive
/// so that the bounding box is [left, right]x[bottom, top].

struct IntBBox {
	/// \brief Default constructor for all zero bounds
	IntBBox(void) : left(0), right(0), bottom(0), top(0), observation_count(0) {}
	/// \brief Count the total number of soundings in the current bounding box
	s32 CountObservations(CountData const*const counts);
	
	u32	left;				///< Minimum column to consider in the box
	u32 right;				///< Maximum column to consider in the box
	u32 bottom;				///< Minimum row to consider in the box
	u32 top;				///< Maimum row to consider in the box
	u32 observation_count;	///< Count of observations in the bounding box
};

/// \brief Output stream operator for a bounding box
std::ostream& operator<<(std::ostream& os, IntBBox const& a);

class DensityTree;	///< Forward declaration of the DensityTree nodes to allow definition of Partition and DensityPartition

/// \class Partition
/// \brief Set of pointers to all the leaf nodes in a particular partition of the area
///
/// A partition is a list of DensityTree nodes that provides the bounding boxes that split up the counts array into as close to
/// equal sounding counts as possible.  This object manages a partition at any stage in the tree that describes all of the possible
/// partitions, and allows the code to accumulate multiple sub-trees so that we can recursively build the list.

class Partition {
public:
	/// \brief Default constructor for an empty partition with zero cost
	Partition(void) : m_aspect_cost(0), m_count_cost(0) {}
	/// \brief Standard constructor to initialise a partition with a single leaf node
	Partition(const DensityTree * const leaf, IntBBox const& box, s32 const nominal_obs);
	
	/// \brief Accumulation operator that concatenates two partial partitions
	Partition& operator+=(Partition const& a);
	
	/// \brief Cost evaluator for the overall cost of the partition currently represented
	inline u32 Cost(void) const { return m_count_cost; }
	
	/// \brief Inspector for the aspect-ratio related cost of the current partition
	inline u32 AspectCost(void) const { return m_aspect_cost; }
	/// \brief Inspector for the sounding-count related cost of the current partition
	inline u32 CountCost(void) const { return m_count_cost; }
	
	/// \typedef const_iterator
	/// \brief Synonym for a const iterator over the list
	typedef std::list<DensityTree const *>::const_iterator const_iterator;
	/// \typedef iterator
	/// \brief Synonym for an iterator over the list
	typedef std::list<DensityTree const *>::iterator iterator;
	
	/// \brief Inspector for  a const iterator to the start of the leaf list
	inline const_iterator begin(void) const { return m_path.begin(); }
	/// \brief Inspector for a const iterator to the end of the leaf list
	inline const_iterator end(void) const { return m_path.end(); }
	/// \brief Inspector for the total number of leaves in the current partition list
	inline u32 size(void) const { return static_cast<u32>( m_path.size() ); }
	
private:
	std::list<DensityTree const *>	m_path;			///< STL list of the leaf nodes currently in the partition
	u32								m_aspect_cost;	///< Cumulative aspect-ratio based cost of the leaves in the current partition
	u32								m_count_cost;	///< Cumulative observation-count based cost of the leaves in the current partition
};

/// \brief Summation operator for partitions, finding the union of the two partial partitions
Partition operator+(Partition const a, Partition const b);

/// \struct CountData;
/// \brief Forward declaration of a computational structure with no public interface
struct CountData;

/// \class DensityTree
/// \brief One node in the partition tree for a given sounding count grid
///
/// This represents a single node in the partition tree for the count grid.  Each node is split into two (not necessarily equal) parts
/// either vertically or horizontally, but in each case we call these 'port' and 'starboard' for convenience; in the vertical case, the
/// port sub-section is the western most, in the horizontal case it is the southern most.  As well as pointers for the sub-divisions,
/// we keep the cost of the best known sub-division to date, and the bounding box for this node.  The constructor calls itself
/// recursively to test all possible sub-divisions of the given bounding box, but rejects all sub-divisions that are provably higher
/// cost that the best known cost provided at construction.  It can also prune out branches of the tree if they are going to cause the
/// total cost (of port and starboard) to exceed the known best cost.  Because of this, the tree is always up to date with the best known
/// partition, although this evolves as the tree is built.

class DensityTree {
public:
	/// \brief Recursive constructor for the optimal partition of the count array
	DensityTree(CountData const* data, u32 const available_cost, u32 const n_sections, s32 const nominal_obs, IntBBox const& bbox);
	/// \brief Default destructor, recursively removing all of the partitions for the tree
	~DensityTree(void);
	
	/// \brief Extract the best available partition in the current tree
	void BestPartition(Partition& p, s32 const nominal_obs) const;
	/// \brief Enumerate all partitions in the current tree
	void EnumeratePartitions(std::vector<Partition>& traces, s32 const nominal_obs);
	/// \brief Convert a given partition into a simple list of bounding boxes that can be passed out to the user
	std::list<IntBBox> *Evaluate(Partition const& p) const;
	/// \brief Inspector for the bounding box represented by this node in the tree
	IntBBox const& BoundingBox(void) const { return m_bbox; }
	/// \brief Inspector for the total number of soundings in the current bounding box
	inline s32 ObservationCount(void) const { return m_bbox.observation_count; }
	/// \brief Generate a representation of the tree on an output stream
	void Dump(std::ostream& os) const;
	/// \brief Inspector for the minimum cost of all possible sub-divisions at this node
	inline u32 MinCost(void) const { return m_lowestCost; }
	
private:
	/// \brief Construct all possible vertical splits of the current node, updating the node with the best one known
	void VerticalSplit(CountData const* data, u32 const split_point, u32 const n_sections, u32 const nominal_obs,
					   std::vector<u32> const& col_partition, u32 const available_cost, IntBBox const& bbox);
	/// \brief Construct all possible horizontal splits of the current node, updating the node with the best one known
	void HorizontalSplit(CountData const* data, u32 const split_point, u32 const n_sections, u32 const nominal_obs,
						 std::vector<u32> const& row_partition, u32 const available_cost, IntBBox const& bbox);
	/// \brief Construct the optimal partition of the bounding box specified into equi-count sections across the width
	void RowPartition(CountData const*const counts, IntBBox const& bbox, u32 const n_sections, std::vector<u32>& partition);
	/// \brief Construct the optimal partition of the bounding box specified into equi-count sections up the height
	void ColPartition(CountData const * const counts, IntBBox const& bbox, u32 const n_sections, std::vector<u32>& partition);
	/// \brief Saturation arithmetic combination of elemental costs from two sub-trees
	u32 CombineCosts(DensityTree const *a, DensityTree const *b);
	
	IntBBox			m_bbox;				///< Bounding box for the current node in the tree
	u32				m_lowestCost;		///< Cost of the best known partition of the tree so far encountered at this node
	DensityTree*	m_portPartition;	///< Port-side sub-tree refinement
	DensityTree*	m_stbdPartition;	///< Starboard side sub-tree refinement
};

/// \class DensityPartition
/// \brief Compute a partition of a given array of counts so that it can be computed in sections
///
/// This provides a helper class to compute a partition of an array of sounding counts into a given
/// number of sections such that each section is as square as possible, and contains as close as possible
/// the same number of soundings.  These conditions assist in making sure that the computations will
/// complete as much as possible at the same time, and that they will minimise any over-spray of lines
/// (and hopefully minimise the number of lines that intersect them).  The solution isn't going to be
/// optimal, but it should assist in getting close.

class DensityPartition {
public:
	/// \brief Default constructor: set up the cumulative counts arrays ready for partitioning
	DensityPartition(u32 const *data, u32 const nCols, u32 const nRows);
	/// \brief Default destructor
	~DensityPartition(void);
	
	/// \brief Compute the optimal partition of a section of the counts array into a given number of sections
	void Repartition(u32 const n_sections, IntBBox const& bbox);
	/// \brief Compute the optimal partitions of the whole counts array into a given number of sections
	void Repartition(u32 const n_sections);
	/// \brief Extract the best known partition of the current tree
	std::list<IntBBox> *BestPartition(u32& cost) const;
	
	// Mostly debugging options, since enumeration of all partitions, or dumping the tree is expensive
	
	/// \brief Find a list of all known partitions of the tree
	void EnumeratePartitions(std::vector<Partition>& traces) const;
	/// \brief Write a text description of the tree to the given output stream
	void DumpTree(std::ostream& os) const;
	/// \brief Inspector for the total number of observations
	u32 ObservationCount(void) const;
	
private:
	CountData				*m_cumulatives;	///< Pointer to the pre-computed cumulative counts arrays
	u32						m_nCols;		///< Width of the counts array
	u32						m_nRows;		///< Height of the counts array
	DensityTree				*m_root;		///< Pointer to the root of the partition tree
	u32						m_nSections;	///< Number of sections in the partition
};

#endif
