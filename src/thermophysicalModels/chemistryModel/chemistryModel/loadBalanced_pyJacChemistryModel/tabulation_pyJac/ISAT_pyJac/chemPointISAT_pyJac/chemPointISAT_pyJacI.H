/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

inline const Foam::chemistryTabulationMethods::ISAT_pyJac&
Foam::chemPointISAT_pyJac::table() const
{
    return table_;
}


inline Foam::label Foam::chemPointISAT_pyJac::nGrowth() const
{
    return nGrowth_;
}


inline Foam::label& Foam::chemPointISAT_pyJac::completeSpaceSize()
{
    return completeSpaceSize_;
}


inline const Foam::label& Foam::chemPointISAT_pyJac::completeSpaceSize() const
{
    return completeSpaceSize_;
}


inline const Foam::scalarField& Foam::chemPointISAT_pyJac::phi() const
{
    return phi_;
}


inline const Foam::scalarField& Foam::chemPointISAT_pyJac::Rphi() const
{
    return Rphi_;
}


inline const Foam::scalarField& Foam::chemPointISAT_pyJac::scaleFactor() const
{
    return scaleFactor_;
}


inline const Foam::scalar& Foam::chemPointISAT_pyJac::tolerance() const
{
    return tolerance_;
}


inline void Foam::chemPointISAT_pyJac::changeTolerance(scalar newTol)
{
    tolerance_ = newTol;
}


inline Foam::binaryNode_pyJac*& Foam::chemPointISAT_pyJac::node()
{
    return node_;
}


inline const Foam::scalarSquareMatrix& Foam::chemPointISAT_pyJac::A() const
{
    return A_;
}


inline Foam::scalarSquareMatrix& Foam::chemPointISAT_pyJac::A()
{
    return A_;
}


inline const Foam::scalarSquareMatrix& Foam::chemPointISAT_pyJac::LT() const
{
    return LT_;
}


inline Foam::scalarSquareMatrix& Foam::chemPointISAT_pyJac::LT()
{
    return LT_;
}


inline Foam::label Foam::chemPointISAT_pyJac::nActive() const
{
    return nActive_;
}




inline void Foam::chemPointISAT_pyJac::increaseNumRetrieve()
{
    this->numRetrieve_++;
}


inline void Foam::chemPointISAT_pyJac::resetNumRetrieve()
{
    this->numRetrieve_ = 0;
}


inline void Foam::chemPointISAT_pyJac::increaseNLifeTime()
{
    this->nLifeTime_++;
}




inline const Foam::label& Foam::chemPointISAT_pyJac::timeTag()
{
    return timeTag_;
}


inline Foam::label& Foam::chemPointISAT_pyJac::lastTimeUsed()
{
    return lastTimeUsed_;
}


inline bool& Foam::chemPointISAT_pyJac::toRemove()
{
    return toRemove_;
}


inline Foam::label& Foam::chemPointISAT_pyJac::maxNumNewDim()
{
    return maxNumNewDim_;
}


inline const Foam::label& Foam::chemPointISAT_pyJac::numRetrieve()
{
    return numRetrieve_;
}


inline const Foam::label& Foam::chemPointISAT_pyJac::nLifeTime()
{
    return nLifeTime_;
}


// ************************************************************************* //
