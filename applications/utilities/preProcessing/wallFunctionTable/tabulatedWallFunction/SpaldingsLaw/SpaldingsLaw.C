/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "SpaldingsLaw.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace tabulatedWallFunctions
    {
        defineTypeNameAndDebug(SpaldingsLaw, 0);
        addToRunTimeSelectionTable
        (
            tabulatedWallFunction,
            SpaldingsLaw,
            dictionary
        );
    }
}

const Foam::label Foam::tabulatedWallFunctions::SpaldingsLaw::maxIters_ = 1000;

const Foam::scalar
    Foam::tabulatedWallFunctions::SpaldingsLaw::tolerance_ = 1e-4;


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::tabulatedWallFunctions::SpaldingsLaw::invertFunction()
{
    // Initialise u+ and R
    scalar Re = 0.0;
    scalar uPlus = 1;

    // Populate the table
    forAll(invertedTable_, i)
    {
        if (invertedTable_.log10())
        {
            Re = pow(10, (i*invertedTable_.dx() + invertedTable_.x0()));
        }
        else
        {
            Re = i*invertedTable_.dx() + invertedTable_.x0();
        }

        // Use latest available u+ estimate
        if (i > 0)
        {
            uPlus = invertedTable_[i-1];
        }

        // Newton iterations to determine u+
        label iter = 0;
        scalar error = great;
        do
        {
            scalar kUPlus = min(kappa_*uPlus, 50);
            scalar A =
                E_*sqr(uPlus)
              + uPlus
               *(exp(kUPlus) - pow3(kUPlus)/6 - 0.5*sqr(kUPlus) - kUPlus - 1);

            scalar f = - Re + A/E_;

            scalar df =
                (
                    2*E_*uPlus
                  + exp(kUPlus)*(kUPlus + 1)
                  - 2.0/3.0*pow3(kUPlus)
                  - 1.5*sqr(kUPlus)
                  - 2*kUPlus
                  - 1
                )/E_;

            scalar uPlusNew = uPlus - f/(df + rootVSmall);
            error = mag((uPlus - uPlusNew)/uPlusNew);
            uPlus = uPlusNew;
        } while (error > tolerance_ && ++iter < maxIters_);

        if (iter == maxIters_)
        {
            WarningInFunction
                << "Newton iterations not converged:" << nl
                << "    iters = " << iter << ", error = " << error << endl;
        }

        // Set new values - constrain u+ >= 0
        invertedTable_[i] = max(0, uPlus);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tabulatedWallFunctions::SpaldingsLaw::SpaldingsLaw
(
    const dictionary& dict,
    const polyMesh& mesh
)
:
    tabulatedWallFunction(dict, mesh, typeName),
    kappa_(readScalar(coeffDict_.lookup("kappa"))),
    E_(readScalar(coeffDict_.lookup("E")))
{
    invertFunction();

    if (debug)
    {
        writeData(Info);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tabulatedWallFunctions::SpaldingsLaw::~SpaldingsLaw()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::tabulatedWallFunctions::SpaldingsLaw::yPlus
(
    const scalar uPlus
) const
{
    scalar kUPlus = min(kappa_*uPlus, 50);

    return
        uPlus
      + 1/E_*(exp(kUPlus) - pow3(kUPlus)/6 - 0.5*sqr(kUPlus) - kUPlus - 1);
}


Foam::scalar Foam::tabulatedWallFunctions::SpaldingsLaw::Re
(
    const scalar uPlus
) const
{
    return uPlus*yPlus(uPlus);
}


void Foam::tabulatedWallFunctions::SpaldingsLaw::writeData(Ostream& os) const
{
    if (invertedTable_.log10())
    {
        os  << "log10(Re), y+, u+:" << endl;
        forAll(invertedTable_, i)
        {
            scalar uPlus = invertedTable_[i];
            scalar Re = ::log10(this->Re(uPlus));
            scalar yPlus = this->yPlus(uPlus);
            os  << Re << ", " << yPlus << ", " << uPlus << endl;
        }
    }
    else
    {
        os  << "Re, y+, u+:" << endl;
        forAll(invertedTable_, i)
        {
            scalar uPlus = invertedTable_[i];
            scalar Re = this->Re(uPlus);
            scalar yPlus = this->yPlus(uPlus);
            os  << Re << ", " << yPlus << ", " << uPlus << endl;
        }
    }
}


// ************************************************************************* //
