/// \file DensityPartition.cpp
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

#include <list>
#include <vector>
#include <limits>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>

#include "stdtypes.h"
#include "DensityPartition.h"

/// \struct CountData
/// \brief Encapsulation of the cumulative counts along rows and down columns
///
/// In order to avoid lots of re-computation, we pre-compute the cumulative sum of counts across each row and
/// down each column in the dataset that we've been given.  We can then compute the sum of observations in any
/// window of the data array by subtracting the cumulative sum of the point prior to the window from that at
/// the end of the window.  In order to make sure that this can always happen, we tack on a zero point at the
/// start of the cumulative sum, making the cumulative arrays one point longer in the summation direction.  The
/// only other mildy interesting thing is that we transpose the sums down the columns so that columns window
/// lookups happen along the direction of the cumulative array (which should give us slightly better cache
/// performance).  This almost-POD is primarily designed to make it easier to pass all of this data around.

struct CountData {
	/// \brief Default constructor for given rectangular array of counts data
	///
	/// Construct the cumulative distributions of observations in rows and columns given the row major (south to north)
	/// source array.
	///
	/// \param data		Pointer to the source data of observation counts per grid cell
	/// \param nCols	Width of the source array in grid cells
	/// \param nRows	Height of the source array in grid cells
	/// \return N/A
	CountData(u32 const *data, u32 const nCols, u32 const nRows)
	: m_origRows(nRows), m_origCols(nCols)
	{
		m_rowCumulatives = new u32[nRows*(nCols+3)];
		m_colCumulatives = new u32[(nRows+3)*nCols];
		
		u32 *p = m_colCumulatives;
		m_nObservations = 0;
		for (u32 c = 0; c < nCols; ++c) {
			*p++ = 0; *p++ = 0;
			for (u32 r = 0; r < nRows; ++r, ++p) {
				*p = *(p-1) + data[c + r*nCols];
			}
			*p = *(p-1);
			m_nObservations += *(p-1);
			++p;
		}
		
		p = m_rowCumulatives;
		for (u32 r = 0; r < nRows; ++r) {
			*p++ = 0; *p++ = 0;
			for (u32 c = 0; c < nCols; ++c, ++p) {
				*p = *(p-1) + data[c + r*nCols];
			}
			*p = *(p-1);
			++p;
		}
		m_rowCumWidth = nCols + 3;
		m_colCumWidth = nRows + 3;
	}
	
	/// \brief Do a single lookup in the cumulative counts for one direction in the array
	///
	/// This is used by either row or column lookup to determine the sum of counts between the start and end point
	/// given.  The count is inclusive of both endpoints.  Note that this assumes that the start and end points are
	/// valid for the array given at construction, and does not attempt to verify that this is the case for each lookup
	/// (essentially for speed).  Note that this only works because the code that computes the cumulative arrays at
	/// construction transposes the column sums so that they are in row-major order; if this changes, the code here would
	/// have to be specialised.  It is the calling routine's job to provide the pointer to the start of the particular
	/// section of the cumulative sum array that is required for the row/column in question (i.e., this code does not
	/// offset from the common base to the right row/column before doing the sum over the window provided).
	///
	/// \param base		Pointer to the cumulative array to use
	/// \param start_pt	First row/column in the window over which to sum
	/// \param end_pt	Last row/column in the window over which to sum
	/// \return Count of observations in the window between the start and end point provided.
	
	inline u32 Lookup(u32 *const base, s32 const start_pt, s32 const end_pt) const
	{
		return base[2+end_pt] - base[2+start_pt-1];
	}
	
	/// \brief Compute the sum of counts in a particular window in a particular row of the array
	///
	/// This computes the sum of all of the counts within a window of a particular row within the cumulative sum array
	/// pre-computed at construction.  The code computes both the sum across the window of columns given, and the sum
	/// over a window one column greater left and right the given window (inner and outer, respectively); the count
	/// over the specified window is returned.  Note that the start and end column of the window are allowed to be
	/// anywhere within the array given, since the pre-computed cumulative arrays are offset to allow for edge effects.
	///	
	/// \param row			Row (in range [0, R-1]) on which to compute
	/// \param start_col	Start column (in range [0, C-1]) at which to start count
	/// \param end_col		End column (in range [0, C-1]) at which to end count
	/// \param inner		Count of all observations in \a row within window [\a start_col, \a end_col]
	/// \param outer		Count of all observations in \a row within window [\a start_col-1, \a end_col+1]
	/// \return	Same as \a inner
	
	inline u32 RowSum(s32 const row, s32 const start_col, s32 const end_col, u32& inner, u32& outer) const
	{
		if (row < 0 || row >= (s32)m_origRows) {
			inner = outer = 0;
			return 0;
		}
		u32 *base = m_rowCumulatives + row*m_rowCumWidth;
		inner = Lookup(base, start_col, end_col);
		outer = Lookup(base, start_col-1, end_col+1);
		return inner;
	}
	
	/// \brief Compute the sum of counts in a particular window in a particular column of the array
	///
	/// This computes the sum of all of the counts within a window of a particular column within the cumulative sum array
	/// pre-computed at construction.  The code computes both the sum across the window of rows given, and the sum
	/// over a window one row greater above and below the given window (inner and outer, respectively); the count
	/// over the specified window is returned.  Note that the start and end row of the window are allowed to be
	/// anywhere within the array given, since the pre-computed cumulative arrays are offset to allow for edge effects.
	///	
	/// \param col			Column (in range [0, C-1]) on which to compute
	/// \param start_row	Start row (in range [0, R-1]) at which to start count
	/// \param end_row		End row (in range [0, R-1]) at which to end count
	/// \param inner		Count of all observations in \a col within window [\a start_row, \a end_row]
	/// \param outer		Count of all observations in \a col within window [\a start_row-1, \a end_row+1]
	/// \return	Same as \a inner
	
	inline u32 ColSum(s32 const col, s32 const start_row, s32 const end_row, u32& inner, u32& outer) const
	{
		if (col < 0 || col >= (s32)m_origCols) {
			inner = outer = 0;
			return 0;
		}
		u32 *base = m_colCumulatives + col*m_colCumWidth;
		inner = Lookup(base, start_row, end_row);
		outer = Lookup(base, start_row-1, end_row+1);
		return inner;
	}
	
	/// Default destructor, simply removing the dynamicaly allocated memory.
	///
	/// \return N/A
	~CountData(void)
	{
		delete m_rowCumulatives;
		delete m_colCumulatives;
	}
	
	u32 *m_rowCumulatives;	///< Row major array of cumulative sums across each row of the counts array
	u32 *m_colCumulatives;	///< Column major array of cumulative sums up each column of the counts array
	u32	m_origRows;			///< Height of the original counts array
	u32	m_origCols;			///< Width of the original counts array
	u32 m_rowCumWidth;		///< Width in columns of the cumulative row array
	u32 m_colCumWidth;		///< Width in rows of the (transposed) cumulative column array
	u32 m_nObservations;	///< Total number of observations found in the counts array
};

/// Output stream operator for an integer-value bounding box.  This generates a decorated output with
/// information about the number of observations in the bounding box, and the number of grid cells that
/// is mostly useful for debugging and human output.
///
/// \param os	Reference for the output stream on which to write the information
/// \param a	Reference for the bounding box to stream
/// \return Reference for the output stream for cascading

std::ostream& operator<<(std::ostream& os, IntBBox const& a)
{
	os << "(" << a.left << ", " << a.right << ") x (" << a.bottom << ", " << a.top << ")";
	if (a.observation_count != 0) {
		u32 cells = (a.right - a.left + 1)*(a.top - a.bottom + 1);
		os << " [" << a.observation_count << " observations, " << cells << " grid cells]";
	}
	return os;
}

/// Output stream for a partition set for the counts array (i.e., a list of the DensityTree (pointer) structures that
/// contain IntBBox objects which are complete and orthogonal over the counts array coverage).  This generates the
/// outputs for all of the bounding boxes sequentially on the output stream using the IntBBox streaming operator, but
/// provides no other decorations.
///
/// \param os	Reference for the output stream on which to write the information
/// \param a	Reference for the partition to stream
/// \return Reference for the output stream for cascading

std::ostream& operator<<(std::ostream& os, Partition const& a)
{
	for (Partition::const_iterator p = a.begin(); p != a.end(); ++p) {
		os << (*p)->BoundingBox() << std::endl;
	}
	return os;
}

/// Count the number of observations in the current bounding box given the cumulative row count for
/// the box.  We could use either row or column cumulatives; using the row cumulative is simply
/// conventional in this case.
///
/// \param counts	Pointer to the pre-computed counts data cache
/// \return Total number of counts in the section of the array indicated by the bounding box limits

s32 IntBBox::CountObservations(CountData const*const counts)
{
	u32 inner, outer;
	observation_count = 0;
	for (u32 row = bottom; row <= top; ++row) {
		observation_count += counts->RowSum(row, left, right, inner, outer);
	}
	return observation_count;
}

/// Cumulative sum operator for partitions.  The cumulation of two partitions simply splices together the lists
/// of DensityTree (pointer) objects, and accumulates the costs associated with the two separate partitions.
///
/// \param a	Reference for the partition to accumulate into the current one
/// \return Reference for the current (and hence updated) partition

Partition& Partition::operator+=(Partition const& a)
{
	m_count_cost += a.m_count_cost;
	m_aspect_cost += a.m_aspect_cost;
	std::list<DensityTree const *> tmp(a.m_path);
	m_path.splice(m_path.end(), tmp);
	return *this;
}

/// Addition operator for two partitions.  This defines the operation in terms of the cumulative operator in
/// the traditional manner.  Note that the routine makes no effort to sort or otherwise arrange the elements of
/// the partition, nor does it attempt to remove duplicates.
///
/// \param a	Reference for the partition to represent first in the output sum
/// \param b	Reference for the partition to represent second in the output sum
/// \return A temporary partition containing the sum of the input partitions

Partition operator+(Partition const a, Partition const b)
{
	Partition c = a;
	return c += b;
}

/// Default constructor for a partition of the counts array.  A partition is a collection of bounding boxes that
/// form a coverage for the input array in the sense that they completely cover the entire array with no overlap
/// between the bounding boxes.  The partition contains pointers to the DensityTree structures that represent the
/// tree nodes in traversal order for these bounding boxes, along with the cumulative cost associated with all of
/// the leaf nodes.  A number of different costs might be associated with a particular partition depending on the
/// code, but the Cost() function can find the appropriate value for use.
///
/// \param leaf					Pointer into the tree node representing the leaf bounding box
/// \param bbox					Bounding box represented in this leaf
/// \param nominal_soundings	Estimate of how many soundings there should ideally be in a leaf node
/// \return N/A

Partition::Partition(const DensityTree * const leaf, IntBBox const& bbox, s32 const nominal_obs)
{
	m_path.push_back(leaf);
	m_aspect_cost = abs(static_cast<s32>(bbox.right - bbox.left + 1) - static_cast<s32>(bbox.top - bbox.bottom + 1));
    m_count_cost = abs(static_cast<s32>(bbox.observation_count - nominal_obs));
}

/// Recursive constructor that builds the entire (pruned) tree from one call on the root node.  The constructor
/// uses three techniques to improve performance: it tries more or less even partitions of the available area first
/// in the hope that a more even distribution, particularly at early stages, will be more likely to be optimal.
/// (Primarily because any errors committed can be amortised over more partitions in the remaining areas.)  Then,
/// it keeps a list of the best current cost of all of the leaves that have been seen so far at each level of the
/// tree and avoids doing the evaluation of any branch of the tree that is not going to be valid (i.e., the part
/// of the branch from the current node that has already been evaluated is over the cost already).  The algorithm
/// also keeps track of the cumulative cost of the current split so that as the split progresses the new branches
/// (at all levels) have to be better and better to evaluate the remainder of the tree.  Finally, the code attempts
/// to improve performance by evaluating any single chunks that it finds (i.e., when the split contains a chunk that
/// has to be a leaf) first.  These are likely to have larger errors overall, and therefore should tell us sooner
/// rather than later whether it's worthwhile to evaluate the remainder of the tree.  The node keeps track only of
/// the best of the current possible splits which strongly limits the amount of memory that is required (this is
/// important since the tree very quickly grows wide and deep as a functon of the number of partition segments
/// required --- so much so that the tree can't fit in any reasonable amount of memory).
///
/// \param data					Cumulative arrays constructed from the counts data array
/// \param available_cost		Maximum cost that this section can accrue and still be plausible
/// \param n_sections			Number of sections into which the current bounding box must be split
/// \param nominal_obs          Expected 'ideal' number of observations in one leaf node of the tree
/// \param bbox					Bounding box for the current section, which is to be split up
/// \return N/A

DensityTree::DensityTree(CountData const* data, u32 const available_cost, u32 const n_sections,
						 s32 const nominal_obs, IntBBox const& bbox)
: m_bbox(bbox)
{
	if (1 == n_sections) {
		// This is the termination of tree construction
		m_bbox.CountObservations(data);
		m_lowestCost = abs(static_cast<s32>(m_bbox.observation_count - nominal_obs));
		m_portPartition = NULL;
		m_stbdPartition = NULL;
		return;
	}
	
	// Since we're going to partition the space, there's no point in computing the nominal cost for
	// this chunk, so we set the lowest cost indicator so that we can estimate it in the loops below
	m_lowestCost = std::numeric_limits<u32>::max();
	m_portPartition = NULL;
	m_stbdPartition = NULL;

	// Process column partitions (i.e., splitting the section vertically)
	std::vector<u32> col_partition;
	ColPartition(data, bbox, n_sections, col_partition);
	std::vector<u32> row_partition;
	RowPartition(data, bbox, n_sections, row_partition);

	s32 section = n_sections/2 - 1, direction;
	if ((n_sections % 2) == 1)
		direction = +1;
	else
		direction = -1;
	
	u32 n = 0;
	u32 target_cost = available_cost;
	do {
		VerticalSplit(data, section, n_sections, nominal_obs, col_partition, target_cost, bbox);
		// If the refinements have resulted in a cost lower than the available cost from above, we do any
		// further refinements with the tighter bound
		if (m_lowestCost < target_cost) target_cost = m_lowestCost;

		HorizontalSplit(data, section, n_sections, nominal_obs, row_partition, target_cost, bbox);
		if (m_lowestCost < target_cost) target_cost = m_lowestCost;
		section += direction*(n+1);
		direction = -direction;
		++n;
	} while (n <= n_sections-2);
}

/// Implement the action to split a given bounding box in the vertical direction, keeping track of the current
/// cost so that we don't evaluate any more of the tree than we need to, keeping track of only the best
/// possible partition, and keeping a running tally of the costs of the partitions.  This calls the constructor
/// for the class recursively to generate all possible combinations of horizontal and vertical splits of each side
/// of the given bounding box, and then evaluates whether each split is better than the best known so far.
///
/// \param data					Cumulative counts array constructed from the counts array
/// \param split_point			Postion within the column partition table at which to split
/// \param n_sections			Number of sections into which the box should be split
/// \param nominal_obs          Ideal number of observations that should appear in a single leaf of the partition
/// \param col_partition		Column partition table: columns at which the most optimal split occurs
/// \param available_cost		Maximum cost available for the split of this bounding box
/// \param bbox					Bounding box to be split
/// \return N/A

void DensityTree::VerticalSplit(CountData const* data, u32 const split_point, u32 const n_sections, u32 const nominal_obs,
								std::vector<u32> const& col_partition, u32 const available_cost, IntBBox const& bbox)
{
	IntBBox port_box = bbox, stbd_box = bbox;
	port_box.right = col_partition[split_point];
	stbd_box.left = col_partition[split_point]+1;
	
	DensityTree *port = NULL;
	DensityTree *stbd = NULL;
	
	if (split_point < n_sections/2) {
		// Since the split-point is to the port-side of the range, we develop the port tree first, and then (if required)
		// the starboard tree.  Since the tree should be shallower on the port-side, we should give ourselves a better chance
		// of not having to develop the starboard-side, which will take longer.
		port = new DensityTree(data, available_cost, split_point+1, nominal_obs, port_box);
		if (port->MinCost() < available_cost) {
			// Looks like we need to evaluate the starboard tree and check whether there is a configuration that would work
			stbd = new DensityTree(data, available_cost - port->MinCost(), n_sections - (split_point + 1), nominal_obs, stbd_box);
		} else {
			// No need to evaluate the starboard tree (since the port tree already exceeds the available cost), and therefore
			// we don't need the port tree either (since it exceeds the available cost).  We can clear up, therefore.
			delete port;
			port = NULL;
		}
	} else {
		// Since the split-point is putting less area on the starboard-side of the split, we develop the starboard side first in
		// the hope that it'll be shallower (and hence faster to evaluate) before seeing whether we need to evaluate the port-side
		// tree (which only happens if there if sufficient cost available after we remove the cost of what we've already done).
		stbd = new DensityTree(data, available_cost, n_sections - (split_point + 1), nominal_obs, stbd_box);
		if (stbd->MinCost() < available_cost) {
			port = new DensityTree(data, available_cost - stbd->MinCost(), split_point + 1, nominal_obs, port_box);
		} else {
			delete stbd;
			stbd = NULL;
		}
	}
	// If we still have a port tree, then we need to check whether there is a partition that we can use: it has to be
	// better than any that we've seen so far to replace the current favourite.
	if (NULL != port) {
		u32 pathCost = CombineCosts(port, stbd);
		if (pathCost < m_lowestCost && pathCost < available_cost) {
			if (NULL != m_portPartition) {
				delete m_portPartition;
				delete m_stbdPartition;
			}
			m_portPartition = port;
			m_stbdPartition = stbd;
			m_lowestCost = pathCost;
		} else {
			delete port;
			delete stbd;
		}
	}	
}

/// Implement the action to split a given bounding box in the horizontal direction, keeping track of the current
/// cost so that we don't evaluate any more of the tree than we need to, keeping track of only the best
/// possible partition, and keeping a running tally of the costs of the partitions.  This calls the constructor
/// for the class recursively to generate all possible combinations of horizontal and vertical splits of each side of
/// the given bounding box, and then evaluates whether each split is better than the best known so far.
///
/// \param data					Cumulative counts array constructed from the counts array
/// \param split_point			Postion within the column partition table at which to split
/// \param n_sections			Number of sections into which the box should be split
/// \param nominal_obs          Ideal number of observations that should appear in a single leaf of the partition
/// \param row_partition		Row partition table: rows at which the most optimal split occurs
/// \param available_cost		Maximum cost available for the split of this bounding box
/// \param bbox					Bounding box to be split
/// \return N/A

void DensityTree::HorizontalSplit(CountData const* data, u32 const split_point, u32 const n_sections, u32 const nominal_obs,
								  std::vector<u32> const& row_partition, u32 const available_cost, IntBBox const& bbox)
{
	IntBBox port_box = bbox, stbd_box = bbox;
	port_box.top = row_partition[split_point];
	stbd_box.bottom = row_partition[split_point]+1;
	
	DensityTree *port = NULL;
	DensityTree *stbd = NULL;
	
	if (split_point < n_sections/2) {
		// Since the split-point is to the port-side of the range, we develop the port tree first, and then (if required)
		// the starboard tree.  Since the tree should be shallower on the port-side (and therefore faster to evaluate), we should
		// give outselves a better chance of not having to develop the starboard-side, which will take longer.
		port = new DensityTree(data, available_cost, split_point+1, nominal_obs, port_box);
		if (port->MinCost() < available_cost) {
			// Looks like we need to evaluate the starboard tree, since there is a partition of the port tree lower
			// than the available cost.
			stbd = new DensityTree(data, available_cost - port->MinCost(), n_sections - (split_point + 1), nominal_obs, stbd_box);
		} else {
			// No need to evaluate the stbd tree (since the port tree already exceed the available cost), and therefore
			// we don't need the port tree either (since it exceeds the available cost).  We can clear up, therefore.
			delete port;
			port = NULL;
		}
	} else {
		// Since the split-point is putting less area on the starboard-side of the split, we develop the starboard side first in
		// the hope that it'll be shallower (and hence faster to evaluate) before seeing whether we need to evaluate the port-side
		// tree (which only happens if there is sufficient cost available after we remove the cost of what we're already done).
		stbd = new DensityTree(data, available_cost, n_sections - (split_point + 1), nominal_obs, stbd_box);
		if (stbd->MinCost() < available_cost) {
			port = new DensityTree(data, available_cost - stbd->MinCost(), split_point + 1, nominal_obs, port_box);
		} else {
			delete stbd;
			stbd = NULL;
		}
	}

	// If we still have a port tree, we need to check if it's lower cost than anything we're seen before to see whether
	// we want to update to it.
	if (NULL != port) {
		u32 pathCost = CombineCosts(port, stbd);
		if (pathCost < m_lowestCost && pathCost < available_cost) {
			if (NULL != m_portPartition) {
				delete m_portPartition;
				delete m_stbdPartition;
			}
			m_portPartition = port;
			m_stbdPartition = stbd;
			m_lowestCost = pathCost;
		} else {
			delete port;
			delete stbd;
		}
	}
}

/// This routine constructs the partition points in the horizontal which would split, as optimally as possible, the
/// specified bounding box into equal quantities of soundings.  The goal is to determine the list of rows that
/// would provide Q, 2Q, 3Q, ..., (N-1)Q soundings to the left of them; we can use the cumulative array to determine
/// all possible splits for the given bounding box by scanning in order.
///
/// \param data			Pointer to the pre-computed count data cache
/// \param bbox			Bounding box to compute over
/// \param n_sections	Number of sections into which the bounding box should be split
/// \param partition	[out] List of all transtion points in the array that give equal distributions
/// \return N/A [partition contains the transition points]

void DensityTree::RowPartition(CountData const*const data, IntBBox const& bbox, u32 const n_sections,
							   std::vector<u32>& partition)
{
	u32 win_rows = 2 + bbox.top - bbox.bottom + 1;
	u32 *counts = new u32[win_rows];
	
	// When we're combining the inner and outer counts, we need to blend; it's easier to compute the blend here
	f64 weight_const = 1.0/M_SQRT2;
	f64 weight_scale = M_SQRT2 - 1.0;
	f64 inner_weight = 1.0 - 1.0/M_SQRT2;
	f64 outer_weight = 1.0/M_SQRT2;
	
	// Generate the count across each row within the window of the bounding box
	u32 total_counts = 0;
	s32 bottom_bound = static_cast<s32>(bbox.bottom), top_bound = static_cast<s32>(bbox.top);
	for (s32 r = bottom_bound-1; r <= top_bound+1; ++r) {
		u32 inner_count, outer_count;
		data->RowSum(r, bbox.left, bbox.right, inner_count, outer_count);
		if (r == bottom_bound-1 || r == top_bound+1) {
			f64 weight = weight_const + weight_scale/(r - (bottom_bound-1) + 2);
			counts[r - (bottom_bound-1)] = static_cast<u32>(weight*(outer_count - inner_count));
		} else {
			counts[r - (bottom_bound-1)] = static_cast<u32>(inner_weight*inner_count + outer_weight*outer_count);
		}
		if (r >= bottom_bound && r <= top_bound) total_counts += inner_count;
	}
	
	// Accumulate the counts across the window, to split into N sections (i.e., with N-1 row values)
	u32 section_target = total_counts / n_sections;
	u32 cumulative_count = counts[0]+counts[1], current_row = 2;
	for (u32 s = 0; s < n_sections-1; ++s) {
		while (current_row < win_rows && (cumulative_count + counts[current_row]/M_SQRT2) < section_target*(s+1))
			cumulative_count += counts[current_row++];
		partition.push_back(current_row - 2 + bottom_bound);
	}
	
	delete[] counts;
}

/// This routine constructs the partition points in the vertical which would split, as optimally as possible, the
/// specified bounding box into equal quantities of soundings.  The goal is to determine the list of columns that
/// would provide Q, 2Q, 3Q, ..., (N-1)Q soundings to the left of them; we can use the cumulative array to determine
/// all possible splits for the given bounding box by scanning in order.
///
/// \param data			Pointer to the pre-computed count data cache
/// \param bbox			Bounding box to compute over
/// \param n_sections	Number of sections into which the bounding box should be split
/// \param partition	[out] List of all transtion points in the array that give equal distributions
/// \return N/A [partition contains the transition points]

void DensityTree::ColPartition(CountData const * const data, IntBBox const& bbox, u32 const n_sections,
							   std::vector<u32>& partition)
{
	u32 win_cols = 2 + bbox.right - bbox.left + 1;
	u32 *counts = new u32[win_cols];
	
	// We need to blend the counts from the interior of the window being split and the exterior surround
	// that we'll also need to process to avoid edge effects.  We compute the blend weights here as much
	// as we can to speed things up later.
	f64 weight_const = 1.0/M_SQRT2;
	f64 weight_scale = M_SQRT2 - 1.0;
	f64 inner_weight = 1.0 - 1.0/M_SQRT2;
	f64 outer_weight = 1.0/M_SQRT2;
	
	// Generate the count down each column within the window of the bounding box
	u32 total_counts = 0;
	s32 left_bound = static_cast<s32>(bbox.left), right_bound = static_cast<s32>(bbox.right);
	for (s32 c = left_bound-1; c <= right_bound+1; ++c) {
		u32 inner_count, outer_count;
		data->ColSum(c, bbox.bottom, bbox.top, inner_count, outer_count);
		if (c == left_bound-1 || c == right_bound+1) {
			f64 weight = weight_const + weight_scale/(c - (left_bound-1) + 2);
			counts[c - (left_bound-1)] = static_cast<u32>(weight*(outer_count - inner_count));
		} else {
			counts[c - (left_bound-1)] = static_cast<u32>(inner_weight*inner_count + outer_weight*outer_count);
		}
		if (c >= left_bound && c <= right_bound) total_counts += inner_count;
	}
	m_bbox.observation_count = total_counts;
	
	// Accumulate the counts across the window, to split into N sections (i.e., with N-1 column values)
	u32 section_target = total_counts / n_sections;
	u32 cumulative_count = counts[0] + counts[1], current_col = 2;
	for (u32 s = 0; s < n_sections-1; ++s) {
		while (current_col < win_cols && (cumulative_count + counts[current_col]/M_SQRT2) < section_target*(s+1))
			cumulative_count += counts[current_col++];
		partition.push_back(current_col - 2 + left_bound);
	}
	
	delete[] counts;
}

/// Compute the combined cost of two branches in the tree (generally representing the split to the left and right
/// of the split point being tested) with saturation arithmetic.  This avoids overflow in U32 due to the 'worst case'
/// cost of the maximum representable value.
///
/// \param a	Pointer to the first sub-tree
/// \param b	Pointer to the second sub-tree
/// \return Combined cost of the two sub-trees, or the maximum possible cost, whichever is smaller

u32 DensityTree::CombineCosts(DensityTree const *a, DensityTree const *b)
{
	u32 overhead_available = std::numeric_limits<u32>::max() - a->m_lowestCost;
	u32 rtn;
	
	if (b->m_lowestCost < overhead_available)
		rtn = a->m_lowestCost + b->m_lowestCost;
	else
		rtn = std::numeric_limits<u32>::max();

	return rtn;
}

/// Default destructor, which just has to remove the port and starboard partition sub-trees (this forms a
/// recursive decent destructor which takes down all of the tree in one call from the root).
///
/// \return N/A
DensityTree::~DensityTree(void)
{
	delete m_portPartition;
	delete m_stbdPartition;
}

/// This walks the tree from the current node and generates a linear list of all of the possible partitions
/// of the area represented in the tree.  This recursively calls itself to refine the port/starboard sub-
/// divisions of the tree.  In fact, since the tree is now critically pruned as we go, there really is only one
/// partition in the tree at any one time, so this is a little superfluous; it was of more use when the tree
/// was fully expanded.
///
/// \param traces       List of the currently known partitions from this point down
/// \param nominal_obs	Ideal number of observations that should occur in a single leaf partition
/// \return N/A

void DensityTree::EnumeratePartitions(std::vector<Partition>& traces, s32 const nominal_obs)
{
	if (NULL == m_portPartition) {
		Partition leaf(this, m_bbox, nominal_obs);
		traces.push_back(leaf);
	} else {
		std::vector<Partition> port_traces;
		m_portPartition->EnumeratePartitions(port_traces, nominal_obs);
		std::vector<Partition> stbd_traces;
		m_stbdPartition->EnumeratePartitions(stbd_traces, nominal_obs);
		
		// Form the Cartessian outer-product of the possible traces to give all possible configurations
		for (u32 t_port = 0; t_port < port_traces.size(); ++t_port) {
			for (u32 t_stbd = 0; t_stbd < stbd_traces.size(); ++t_stbd) {
				Partition combo = port_traces[t_port] + stbd_traces[t_stbd];
				traces.push_back(combo);
			}
		}
	}
}

/// Find the best available partition in the curent tree from this point.  Note that because the tree is
/// critically pruned, the only partition in the tree is the best one, so this basically walks the tree
/// and accumulates the results.
///
/// \param p			[out] Space to report the optimal partition
/// \param nominal_obs	Ideal number of observations that should occur in a single leaf partition
/// \return N/A

void DensityTree::BestPartition(Partition& p, s32 const nominal_obs) const
{
	if (NULL == m_portPartition) {
		// We're a leaf node, so we need to add ourselves to the partition
		Partition leaf(this, m_bbox, nominal_obs);
		p += leaf;
	} else {
		// We're somewhere in the hierarchy, and we need to propagate down to the leaves, using the best path
		// only
		m_portPartition->BestPartition(p, nominal_obs);
		m_stbdPartition->BestPartition(p, nominal_obs);
	}
}

/// Convert the partiton from internal format to a simple list of IntBBox structures.  This allows it to be
/// sent to users in more readily usable fashion.
///
/// \param p	Reference for the partition to convert
/// \return Pointer to the simple STL list of the IntBBox objects converted from internal format

std::list<IntBBox> *DensityTree::Evaluate(Partition const& p) const
{
	std::list<IntBBox> *rtn = new std::list<IntBBox>;
	
	for (Partition::const_iterator leaf = p.begin(); leaf != p.end(); ++leaf) {
		rtn->push_back((*leaf)->m_bbox);
	}
	return rtn;
}

/// Dump the entire tree in order from the current node to the leaves, printing the bounding box information as we go.
/// This is primarily useful for debug.
///
/// \param os	Reference for the output stream into which to dump the data
/// \return N/A

void DensityTree::Dump(std::ostream& os) const
{
	os << "Node " << this << " with BBox " << m_bbox << std::endl;
	if (NULL != m_portPartition) {
		os << "Children: (" << m_portPartition << ", " << m_stbdPartition << ")" << std::endl;
		m_portPartition->Dump(os);
		m_stbdPartition->Dump(os);
	} else {
		os << "\tTERMINAL NODE" << std::endl;
	}
}

/// Default constructor for the density partition.  Since we don't know, to start with, what the number of
/// partitions is going to be, we can't do anything more at this stage that pre-compute the cumulative counts
/// arrays.
///
/// \param data		Pointer to the row-major counts data array
/// \param nCols	Width of the counts array
/// \param nRows	Height of the counts array
/// \return N/A

DensityPartition::DensityPartition(u32 const *data, u32 const nCols, u32 const nRows)
: m_nCols(nCols), m_nRows(nRows), m_root(NULL), m_nSections(0)
{
	m_cumulatives = new CountData(data, m_nCols, m_nRows);
}

/// Default destructor: remove the cached cumulative counts array, and then delete the partition tree.
///
/// \return N/A

DensityPartition::~DensityPartition(void)
{
	delete m_cumulatives;
	delete m_root;
}

/// Construct the optimal partition of the configured count array for a certain number of sections.  This is
/// simply a wrapper around the routine that includes bounding box information which defaults the bounding box
/// to the entire size of the input array.
///
/// \param n_sections	Number of sections into which to break the counts array
/// \return N/A

void DensityPartition::Repartition(u32 const n_sections)
{
	IntBBox bbox;
	bbox.left = 0; bbox.bottom = 0;
	bbox.right = m_nCols - 1; bbox.top = m_nRows - 1;
	Repartition(n_sections, bbox);
}

/// Construct the optimal partition of a specific sub-section of the counts array into the given number of sections.
/// This simply constructs the partition tree root node with appropriate limits to find the best partition.
///
/// \param n_sections	Number of sections into which the bounding box should be split
/// \param bbox			Region of the original array to split up
/// \return N/A

void DensityPartition::Repartition(u32 const n_sections, IntBBox const& bbox)
{
	m_nSections = n_sections;
	if (NULL != m_root) delete m_root;
	s32 nominal_obs = m_cumulatives->m_nObservations / n_sections;
	m_root = new DensityTree(m_cumulatives, std::numeric_limits<u32>::max(), n_sections, nominal_obs, bbox);
}

/// Find the optimal partition from all those stored in the tree of possibilities.  Because the tree is critically
/// pruned, this resolves to simply walking the only available partition.  This routine essentially just does the
/// edge computations to work out the nominal number of observations per leaf, and convert the tree walk to just the
/// bounding box information that can be passed back out to the users for general use.
///
/// \param cost	[out] Cost of the partition being returned
/// \return Pointer to the STL list of the bounding boxes for the partition in traversal order

std::list<IntBBox> *DensityPartition::BestPartition(u32& cost) const
{
	Partition p;
	m_root->BestPartition(p, m_cumulatives->m_nObservations/m_nSections);
	std::list<IntBBox> *rtn = m_root->Evaluate(p);
	cost = p.Cost();

	return rtn;
}

/// Find all possible partitions represented in the current tree.  This is essentially a wrapper call for the
/// tree's enumeration call, but since the tree is critically pruned, there really is only one parition that's
/// available.  This routine was more useful before the tree pruning was implemented.
///
/// \param traces	[out] Referece for the STL vector in which to store all possible traces
/// \return N/A

void DensityPartition::EnumeratePartitions(std::vector<Partition>& traces) const
{
	m_root->EnumeratePartitions(traces, m_cumulatives->m_nObservations/m_nSections);
}

/// Write a representation of the tree to the output stream provided.  This writes as much information as possible
/// about the internal structure of the tree, which can be used to reconstruct the internal structure if required.
/// This is useful primarily for debugging.
///
/// \param os	Reference for the output stream on which to generate the tree information
/// \return N/A

void DensityPartition::DumpTree(std::ostream& os) const
{
	os << "--------------------------" << std::endl;
	os << "Root Node: " << m_root << std::endl;
	os << "Total Number of Sections: " << m_nSections << std::endl;
	m_root->Dump(os);
}

/// Inspector for the total number of soundings that were seen in the counts array
///
/// \return Count of the soundings in the input array

u32 DensityPartition::ObservationCount(void) const
{
	return m_cumulatives->m_nObservations;
}
