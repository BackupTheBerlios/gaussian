/* WARANTY NOTICE AND COPYRIGHT
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

Copyright (C) Michael J. Meyer

matmjm@mindspring.com
spyqqqdia@yahoo.com

*/

#ifndef gpr_typedefsmacros_h
#define gpr_typedefsmacros_h

#include <cassert>
#include <iostream>                         // once and for all
using std::cout;
using std::cerr;
using std::cin;
using std::endl;


#define GPR_BEGIN_NAMESPACE(name) namespace name {
#define GPR_END_NAMESPACE(name) }



#define PI 3.14159265358979323846
#define TWO_PI 6.28318531
#define LNPI 1.14472988584940016
#define LNSQRT2PI 0.9189385332046727
#define PNT68 0.6796875
#define SQRT_TWO 1.414213562373095049

#define GR 0.6180339887    // golden ratio, 1/GR=1+GR
#define SMALL 30           // small matrix, uses cache for multiplication
#define SIGMA 0.01         // what we add to the diagonal of the kernel matrix
                           // (K(s_i,s_j)) to ensure positive definiteness

// the basic scalar type
typedef double Real;
// Real function of one real variable
typedef Real (*RealFunction)(Real);




#endif
