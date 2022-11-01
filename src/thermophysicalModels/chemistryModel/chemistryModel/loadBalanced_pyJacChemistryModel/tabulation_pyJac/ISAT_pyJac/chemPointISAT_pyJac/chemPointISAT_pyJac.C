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

#include "chemPointISAT_pyJac.H"
#include "ISAT_pyJac.H"
#include "odeChemistryModel.H"
#include "SVD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Defined as static to be able to dynamically change it during simulations
// (all chemPoints refer to the same object)
Foam::scalar Foam::chemPointISAT_pyJac::tolerance_;

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::chemPointISAT_pyJac::qrDecompose
(
    const label nCols,
    scalarSquareMatrix& R
)
{
    scalarField c(nCols);
    scalarField d(nCols);
    scalar scale, sigma, sum;

    for (label k=0; k<nCols-1; k++)
    {
        scale = 0;
        for (label i=k; i<nCols; i++)
        {
            scale=max(scale, mag(R(i, k)));
        }
        if (scale == 0)
        {
            c[k] = d[k] = 0;
        }
        else
        {
            for (label i=k; i<nCols; i++)
            {
                R(i, k) /= scale;
            }
            sum = 0;
            for (label i=k; i<nCols; i++)
            {
                sum += sqr(R(i, k));
            }
            sigma = sign(R(k, k))*sqrt(sum);
            R(k, k) += sigma;
            c[k] = sigma*R(k, k);
            d[k] = -scale*sigma;
            for (label j=k+1; j<nCols; j++)
            {
                sum = 0;
                for ( label i=k; i<nCols; i++)
                {
                    sum += R(i, k)*R(i, j);
                }
                scalar tau = sum/c[k];
                for ( label i=k; i<nCols; i++)
                {
                    R(i, j) -= tau*R(i, k);
                }
            }
        }
    }

    d[nCols-1] = R(nCols-1, nCols-1);

    // form R
    for (label i=0; i<nCols; i++)
    {
        R(i, i) = d[i];
        for ( label j=0; j<i; j++)
        {
            R(i, j)=0;
        }
    }
}


void Foam::chemPointISAT_pyJac::qrUpdate
(
    scalarSquareMatrix& R,
    const label n,
    const scalarField& u,
    const scalarField& v
)
{
    label k;

    scalarField w(u);
    for (k=n-1; k>=0; k--)
    {
        if (w[k] != 0)
        {
            break;
        }
    }

    if (k < 0)
    {
        k = 0;
    }

    for (label i=k-1; i>=0; i--)
    {
        rotate(R, i, w[i], -w[i+1], n);
        if (w[i] == 0)
        {
            w[i] = mag(w[i+1]);
        }
        else if (mag(w[i]) > mag(w[i+1]))
        {
            w[i] = mag(w[i])*sqrt(1.0 + sqr(w[i+1]/w[i]));
        }
        else
        {
            w[i] = mag(w[i+1])*sqrt(1.0 + sqr(w[i]/w[i+1]));
        }
    }

    for (label i=0; i<n; i++)
    {
        R(0, i) += w[0]*v[i];
    }

    for (label i=0; i<k; i++)
    {
        rotate(R, i, R(i, i), -R(i+1, i), n);
    }
}


void Foam::chemPointISAT_pyJac::rotate
(
    scalarSquareMatrix& R,
    const label i,
    const scalar a,
    const scalar b,
    label n
)
{
    scalar c, fact, s, w, y;

    if (a == 0)
    {
        c = 0;
        s = (b >= 0 ? 1 : -1);
    }
    else if (mag(a) > mag(b))
    {
        fact = b/a;
        c = sign(a)/sqrt(1.0 + sqr(fact));
        s = fact*c;
    }
    else
    {
        fact = a/b;
        s = sign(b)/sqrt(1.0 + sqr(fact));
        c = fact*s;
    }

    for (label j=i;j<n;j++)
    {
        y = R(i, j);
        w = R(i+1, j);
        R(i, j) = c*y-s*w;
        R(i+1, j) = s*y+c*w;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::chemPointISAT_pyJac::chemPointISAT_pyJac
(
    chemistryTabulationMethods::ISAT_pyJac& table,
    const scalarField& phi,
    const scalarField& Rphi,
    const scalarSquareMatrix& A,
    const scalarField& scaleFactor,
    const scalar tolerance,
    const label completeSpaceSize,
    const label nActive,
    const dictionary& coeffsDict,
    binaryNode_pyJac* node
)
:
    table_(table),
    phi_(phi),
    Rphi_(Rphi),
    A_(A),
    scaleFactor_(scaleFactor),
    node_(node),
    completeSpaceSize_(completeSpaceSize),
    nGrowth_(0),
    nActive_(nActive),
    timeTag_(table.timeSteps()),
    lastTimeUsed_(table.timeSteps()),
    toRemove_(false),
    maxNumNewDim_(coeffsDict.lookupOrDefault("maxNumNewDim",0)),
    printProportion_(coeffsDict.lookupOrDefault("printProportion",false)),
    numRetrieve_(0),
    nLifeTime_(0)
{
    tolerance_ = tolerance;

    iddeltaT_ = completeSpaceSize - 1;
    scaleFactor_[iddeltaT_] *= phi_[iddeltaT_]/tolerance_;

    idT_ = completeSpaceSize - 3;
    idp_ = completeSpaceSize - 2;

    const label reduOrCompDim = completeSpaceSize;
    // SVD decomposition A = U*D*V^T
    SVD svdA(A);

    scalarDiagonalMatrix D(reduOrCompDim);
    const scalarDiagonalMatrix& S = svdA.S();

    // Replace the value of vector D by max(D, 1/2), first ISAT paper
    for (label i=0; i<reduOrCompDim; i++)
    {
        D[i] = max(S[i], 0.5);
    }

    // Rebuild A with max length, tol and scale factor before QR decomposition
    scalarRectangularMatrix Atilde(reduOrCompDim);

    // Result stored in Atilde
    multiply(Atilde, svdA.U(), D, svdA.V().T());

    for (label i=0; i<reduOrCompDim; i++)
    {
        for (label j=0; j<reduOrCompDim; j++)
        {
            label compi = i;
            // SF*A/tolerance
            // (where SF is diagonal with inverse of scale factors)
            // SF*A is the same as dividing each line by the scale factor
            // corresponding to the species of this line
            Atilde(i, j) /= (tolerance*scaleFactor[compi]);
        }
    }

    // The object LT_ (the transpose of the Q) describe the EOA, since we have
    // A^T B^T B A that should be factorised into L Q^T Q L^T and is set in the
    // qrDecompose function
    LT_ = scalarSquareMatrix(Atilde);

    qrDecompose(reduOrCompDim, LT_);
}


Foam::chemPointISAT_pyJac::chemPointISAT_pyJac
(
    Foam::chemPointISAT_pyJac& p
)
:
    table_(p.table_),
    phi_(p.phi()),
    Rphi_(p.Rphi()),
    LT_(p.LT()),
    A_(p.A()),
    scaleFactor_(p.scaleFactor()),
    node_(p.node()),
    completeSpaceSize_(p.completeSpaceSize()),
    nGrowth_(p.nGrowth()),
    nActive_(p.nActive()),
    timeTag_(p.timeTag()),
    lastTimeUsed_(p.lastTimeUsed()),
    toRemove_(p.toRemove()),
    maxNumNewDim_(p.maxNumNewDim()),
    numRetrieve_(0),
    nLifeTime_(0)
{
    tolerance_ = p.tolerance();

    idT_ = completeSpaceSize() - 3;
    idp_ = completeSpaceSize() - 2;
    iddeltaT_ = completeSpaceSize() - 1;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::chemPointISAT_pyJac::inEOA(const scalarField& phiq)
{
    const scalarField dphi(phiq - phi());
    const label dim = completeSpaceSize() - 3;

    scalar epsTemp = 0;
    List<scalar> propEps(completeSpaceSize(), scalar(0));

    for (label i=0; i<completeSpaceSize()-3; i++)
    {
        scalar temp = 0;

        // When mechanism reduction is inactive OR on active species multiply L
        // by dphi to get the distance in the active species direction else (for
        // inactive species), just multiply the diagonal element and dphi
        if(true)
        {
            const label si = i;

            for (label j=si; j<dim; j++)// LT is upper triangular
            {
                const label sj = j;

                temp += LT_(si, j)*dphi[sj];
            }

            temp += LT_(si, dim)*dphi[idT_];
            temp += LT_(si, dim+1)*dphi[idp_];
            temp += LT_(si, dim+2)*dphi[iddeltaT_];
        }
        else
        {
            temp = dphi[i]/(tolerance_*scaleFactor_[i]);
        }

        epsTemp += sqr(temp);

        if (printProportion_)
        {
            propEps[i] = temp;
        }
    }
    // Temperature
    epsTemp +=
        sqr
        (
            LT_(dim, dim)*dphi[idT_]
           +LT_(dim, dim+1)*dphi[idp_]
           +LT_(dim, dim+2)*dphi[iddeltaT_]
        );

    // Pressure
    epsTemp +=
        sqr
        (
            LT_(dim+1, dim+1)*dphi[idp_]
           +LT_(dim+1, dim+2)*dphi[iddeltaT_]
        );

    epsTemp += sqr(LT_[dim+2][dim+2]*dphi[iddeltaT_]);

    if (printProportion_)
    {
        propEps[idT_] = sqr
        (
            LT_(dim, dim)*dphi[idT_]
          + LT_(dim, dim+1)*dphi[idp_]
        );

        propEps[idp_] =
            sqr(LT_(dim+1, dim+1)*dphi[idp_]);

        propEps[iddeltaT_] =
            sqr(LT_[dim+2][dim+2]*dphi[iddeltaT_]);
    }

    if (sqrt(epsTemp) > 1 + tolerance_)
    {
        if (printProportion_)
        {
            scalar max = -1;
            label maxIndex = -1;
            for (label i=0; i<completeSpaceSize(); i++)
            {
                if(max < propEps[i])
                {
                    max = propEps[i];
                    maxIndex = i;
                }
            }
            word propName;
            if (maxIndex >= completeSpaceSize() - 3)
            {
                if (maxIndex == idT_)
                {
                    propName = "T";
                }
                else if (maxIndex == idp_)
                {
                    propName = "p";
                }
                else if (maxIndex == iddeltaT_)
                {
                    propName = "deltaT";
                }
            }
            else
            {
                propName = table_.chemistry().Y()[maxIndex].member();
            }

            Info<< "Direction maximum impact to error in ellipsoid: "
                << propName << endl;

            Info<< "Proportion to the total error on the retrieve: "
                << max/(epsTemp+small) << endl;
        }
        return false;
    }
    else
    {
        return true;
    }
}


bool Foam::chemPointISAT_pyJac::checkSolution
(
    const scalarField& phiq,
    const scalarField& Rphiq
)
{
    scalar eps2 = 0;
    const scalarField dR(Rphiq - Rphi());
    const scalarField dphi(phiq - phi());
    const scalarField& scaleFactorV(scaleFactor());
    const scalarSquareMatrix& Avar(A());
    scalar dRl = 0;

    // Since we build only the solution for the species, T and p are not
    // included
    for (label i=0; i<completeSpaceSize()-3; i++)
    {
        dRl = 0;
     
        for (label j=0; j<completeSpaceSize(); j++)
        {
            dRl += Avar(i, j)*dphi[j];
        }
    
        eps2 += sqr((dR[i]-dRl)/scaleFactorV[i]);
    }

    eps2 = sqrt(eps2);
    if (eps2 > tolerance())
    {
        return false;
    }
    else
    {
        // if the solution is in the ellipsoid of accuracy
        return true;
    }
}


bool Foam::chemPointISAT_pyJac::grow(const scalarField& phiq)
{
    const scalarField dphi(phiq - phi());
    const label dim = completeSpaceSize();

    // beginning of grow algorithm
    scalarField phiTilde(dim, 0);
    scalar normPhiTilde = 0;
    // p' = L^T.(p-phi)

    for (label i=0; i<dim; i++)
    {
        for (label j=i; j<dim-3; j++)// LT is upper triangular
        {
            const label sj = j;

            phiTilde[i] += LT_(i, j)*dphi[sj];
        }

        phiTilde[i] += LT_(i, dim-3)*dphi[idT_];
        phiTilde[i] += LT_(i, dim-3+1)*dphi[idp_];
        phiTilde[i] += LT_(i, dim-3+2)*dphi[iddeltaT_];

        normPhiTilde += sqr(phiTilde[i]);
    }

    const scalar invSqrNormPhiTilde = 1.0/normPhiTilde;
    normPhiTilde = sqrt(normPhiTilde);

    // gamma = (1/|p'| - 1)/|p'|^2
    const scalar gamma = (1/normPhiTilde - 1)*invSqrNormPhiTilde;
    const scalarField u(gamma*phiTilde);
    scalarField v(dim, 0);

    for (label i=0; i<dim; i++)
    {
        for (label j=0; j<=i;j++)
        {
            v[i] += phiTilde[j]*LT_(j, i);
        }
    }

    qrUpdate(LT_,dim, u, v);
    nGrowth_++;

    return true;
}


// ************************************************************************* //