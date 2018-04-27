/*! \file chrt_stdtypes.h
 *  \brief Standard types to allow for ready swap from machine to machine
 *
 * This file is used to provide standard type names that indicate actual size, rather than
 * rely on whatever the system gives you.
 */

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

#ifndef __CHRT_STDTYPES_H__
#define __CHRT_STDTYPES_H__

#ifdef __cplusplus
extern "C" {
#endif

typedef unsigned char		u8,		/*!< 8-bit unsigned integer */
							*u8p,	/*!< Pointer to 8-bit unsigned integer */
							byte;	/*!< Synonym for u8 */
typedef signed char			s8,		/*!< 8-bit signed integer */
							*s8p;	/*!< Pointer to 8-bit signed integer */
typedef unsigned short		u16,	/*!< 16-bit unsigned integer */
							*u16p;	/*!< Pointer to 16-bit unsigned integer */
typedef signed short		s16,	/*!< 16-bit signed integer */
							*s16p;	/*!< Pointer to 16-bit signed integer */
typedef unsigned int		u32,	/*!< 32-bit unsigned integer */
							*u32p;	/*!< Pointer to 32-bit unsigned integer */
typedef signed int			s32,	/*!< 32-bit signed integer */
							*s32p;	/*!< Pointer to 32-bit signed integer */
typedef unsigned long long	u64,	/*!< 64-bit unsigned integer */
							*u64p;	/*!< Pointer to 64-bit unsigned integer */
typedef signed long long	s64,	/*!< 64-bit signed integer */
							*s64p;	/*!< Pointer to 64-bit signed integer */
typedef float				f32,	/*!< Single-precision (32-bit) IEEE floating point number */
							*f32p;	/*!< Pointer to single-precision (32-bit) IEEE floating point number */
typedef double				f64,	/*!< Double-precision (64-bit) IEEE floating point number */
							*f64p;	/*!< Pointer to double-precision (64-bit) IEEE floating point number */

#ifdef __cplusplus
}
#endif

#endif	/* __CHRT_STDTYPES_H__ */
