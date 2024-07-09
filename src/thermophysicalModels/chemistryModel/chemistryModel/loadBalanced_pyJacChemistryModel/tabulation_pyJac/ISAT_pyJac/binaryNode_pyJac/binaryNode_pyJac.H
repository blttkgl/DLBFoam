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

Class
    Foam::binaryNode_pyJac

Description
    Node of the binary tree

SourceFile
    binaryNode_pyJac.C

\*---------------------------------------------------------------------------*/

#ifndef binaryNode_pyJac_H
#define binaryNode_pyJac_H

#include "chemPointISAT_pyJac.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class binaryNode_pyJac Declaration
\*---------------------------------------------------------------------------*/

class binaryNode_pyJac
{

public:

    //- Element on the left
    chemPointISAT_pyJac* leafLeft_;

    //- Element on the right
    chemPointISAT_pyJac* leafRight_;

    //- Node which follows on the left
    binaryNode_pyJac* nodeLeft_;

    //- Node which follows on the right
    binaryNode_pyJac* nodeRight_;

    //- Parent node
    binaryNode_pyJac* parent_;

    scalarField v_;

    scalar a_;

    //- Compute vector v:
    //  Let E be the ellipsoid which covers the region of accuracy of
    //  the left leaf (previously defined). E is described by
    //  E={phi| ||L^T.(phi-phi0)|| <= 1}, (see chemPoint for more details).
    //  let E' be the transformation of E in a space where E' is a hypersphere
    //  centered at the origin, in this space y=L^T.(phi-phi0) and then
    //  E'={y| ||y||<=1}
    //  let u be the unit vector joining the center of E' and the newly added
    //  composition point in the transformed space
    //  (y2=L^T.(phiq-phi0)),u = y2/||y2|
    //  Then the hyperplane separating the two points is defined as the
    //  perpendicular bisector of the segment linking 0 to y2
    //      H' = {y| u^T.(y-yh) = 0},
    //  where yh = y2/2.
    //  In the original composition space, the hyperplane H is defined by
    //      H = {y| v^T(phi-phih) = 0},
    //  where phih = phi0 + L^-T.yh = (phi0 + phiq) / 2 and v is
    //            L.L^T (phiq-phi0)
    //      v = -------------------- .
    //          ||L.L^T (phiq-phi0)||
    //
    //  Note :  v is not normalised in this implementation since it is used
    //              on both side of an equality to know if one should go on the
    //              left or the right in the binary tree
    //  Parameter:
    //      elementLeft : chemPoint of the left element
    //      elementRight: chemPoint of the right element
    //       v : empty scalar field to store v
    //  Returnq: void (v is stored in the empty scalarField)
    void calcV
    (
        const chemPointISAT_pyJac& elementLeft,
        const chemPointISAT_pyJac& elementRight,
        scalarField& v
    );

    //- Compute a the product v^T.phih, with phih = (phi0 + phiq)/2.
    //  When travelling in the binary tree,
    //  to know in which part of the composition space the query point 'phi'
    //  belongs to, v^T.phi is performed. If the result is "> a" then it belongs
    //  to the right part (where phiq is), otherwise it belongs to the left
    // part (where phi0 is).
    inline scalar calcA
    (
        const chemPointISAT_pyJac& elementLeft,
        const chemPointISAT_pyJac& elementRight
    );

    // Constructors

        //- Construct null
        binaryNode_pyJac();

        //- Construct from components
        binaryNode_pyJac
        (
            chemPointISAT_pyJac* elementLeft,
            chemPointISAT_pyJac* elementRight,
            binaryNode_pyJac* parent
        );


    // Member Functions

        //- Access

            inline chemPointISAT_pyJac*& leafLeft()
            {
                return leafLeft_;
            }

            inline chemPointISAT_pyJac*& leafRight()
            {
                return leafRight_;
            }

            inline binaryNode_pyJac*& nodeLeft()
            {
                return nodeLeft_;
            }

            inline binaryNode_pyJac*& nodeRight()
            {
                return nodeRight_;
            }

            inline binaryNode_pyJac*& parent()
            {
                return parent_;
            }

        //- Topology

            inline const scalarField& v() const
            {
                return v_;
            }

            inline const scalar& a() const
            {
                return a_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "binaryNode_pyJacI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
