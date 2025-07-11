/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "Reaction2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*
template<class Reaction2Thermo>
inline word& Reaction2<Reaction2Thermo>::name() // can't return a reference of a const member with a non-const-return-valued-function
{
    return name_;
}
*/

template<class Reaction2Thermo>
inline const word& Reaction2<Reaction2Thermo>::name() const
{
    return name_;
}


template<class Reaction2Thermo>
inline const List<typename Reaction2<Reaction2Thermo>::specieCoeffs>&
Reaction2<Reaction2Thermo>::lhs() const
{
    return lhs_;
}


template<class Reaction2Thermo>
inline const List<typename Reaction2<Reaction2Thermo>::specieCoeffs>&
Reaction2<Reaction2Thermo>::rhs() const
{
    return rhs_;
}


template<class Reaction2Thermo>
inline const ControllingTemperatureType&
Reaction2<Reaction2Thermo>::controlT() const
{
    return controlT_;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Reaction2Thermo>
inline Ostream& operator<<(Ostream& os, const Reaction2<Reaction2Thermo>& r)
{
    OStringStream reaction;
    os << r.reactionStr(reaction)<< token::END_STATEMENT <<nl;
   return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
