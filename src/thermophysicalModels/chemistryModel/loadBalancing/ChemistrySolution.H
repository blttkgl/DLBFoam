/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | DLBFoam: Dynamic Load Balancing
   \\    /   O peration     | for fast reactive simulations
    \\  /    A nd           |
     \\/     M anipulation  | 2020, Aalto University, Finland
-------------------------------------------------------------------------------
License
    This file is part of DLBFoam library, derived from OpenFOAM.

    https://github.com/Aalto-CFD/DLBFoam

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

Class
    Foam::ChemistrySolution

Description
    A small object containing everything required for updating the reaction rate
    and the chemistry time step. These are passed around in the load balancer.

\*---------------------------------------------------------------------------*/

#ifndef ChemistrySolution_H
#define ChemistrySolution_H

#include "volFields.H"

namespace Foam
{

struct ChemistrySolution
{

    ChemistrySolution() = default;

    ChemistrySolution(label nspecie)
        : rr(nspecie, 0.0), deltaTChem(0.0), cellid(0), rhoi(0.0), retrieved(false)
    {
    }

    bool operator==(const ChemistrySolution& rhs) const
    {
        return false;
    }

    bool operator!=(const ChemistrySolution& rhs) const
    {
        return !(*this == rhs);
    }

    scalarField rr; // =   (Y_{i+1} - Y_{i}) * rhoi / deltaT;
    scalar deltaTChem;
    scalar cpuTime;
    label cellid;
    scalar rhoi;
    bool retrieved;
};

//- Serialization for send
static inline Ostream& operator<<(Ostream& os, const ChemistrySolution& s)
{
    os << s.rr;
    os << s.deltaTChem;
    os << s.cpuTime;
    os << s.cellid;
    os << s.rhoi;
    os << s.retrieved;
    return os;
}

//- Get a serialized solution from IStream
static inline Istream& operator>>(Istream& is, ChemistrySolution& s)
{
    is >> s.rr;
    is >> s.deltaTChem;
    is >> s.cpuTime;
    is >> s.cellid;
    is >> s.rhoi;
    is >> s.retrieved;
    return is;
}

} // namespace Foam

#endif
