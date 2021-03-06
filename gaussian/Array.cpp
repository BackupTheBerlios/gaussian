/***************************************************************************
                          Array.cpp  -  description
                             -------------------
    begin                : Sat Jun 19 2004
    copyright            : (C) 2004 by 
    email                : 
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "Array.h"

GPR_BEGIN_NAMESPACE(Gaussian)

ostream&
operator << (ostream& os, const RealArray1D& y)
{ y.printSelf(os); }

GPR_END_NAMESPACE(Gaussian)