/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "general.H"
#include "addToRunTimeSelectionTable.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace tabulatedWallFunctions
    {
        defineTypeNameAndDebug(general, 0);
        addToRunTimeSelectionTable
        (
            tabulatedWallFunction,
            general,
            dictionary
        );
    }

    template<>
    const char* Foam::NamedEnum
    <
        Foam::tabulatedWallFunctions::general::interpolationType,
        1
    >::names[] =
    {
        "linear"
    };

}

const
Foam::NamedEnum<Foam::tabulatedWallFunctions::general::interpolationType, 1>
    Foam::tabulatedWallFunctions::general::interpolationTypeNames_;


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::tabulatedWallFunctions::general::invertTable()
{
    scalarList Rey(uPlus_.size(), 0.0);

    // Calculate Reynolds number
    forAll(uPlus_, i)
    {
        Rey[i] = yPlus_[i]*uPlus_[i];
        if (invertedTable_.log10())
        {
            Rey[i] = ::log10(max(ROOTVSMALL, Rey[i]));
        }
    }

    // Populate the U+ vs Re table
    forAll(invertedTable_, i)
    {
        scalar Re = i*invertedTable_.dx() + invertedTable_.x0();
        invertedTable_[i] = max(0, interpolate(Re, Rey, uPlus_));
    }
}


Foam::scalar Foam::tabulatedWallFunctions::general::interpolate
(
    const scalar xi,
    const scalarList& x,
    const scalarList& fx
) const
{
    switch (interpType_)
    {
        case itLinear:
        {
            if (xi <= x[0])
            {
                return fx[0];
            }
            else if (xi >= x.last())
            {
                return fx.last();
            }
            else
            {
                label i2 = 0;
                while (x[i2] < xi)
                {
                    i2++;
                }
                label i1 = i2 - 1;

                return (xi - x[i1])/(x[i2] - x[i1])*(fx[i2] - fx[i1]) + fx[i1];
            }

            break;
        }
        default:
        {
            FatalErrorIn
            (
                "tabulatedWallFunctions::general::interpolate"
                "("
                    "const scalar, "
                    "const scalarList&, "
                    "const scalarList&"
                ")"
            )   << "Unknown interpolation method" << nl
                << abort(FatalError);
        }
    }

    return 0.0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tabulatedWallFunctions::general::general
(
    const dictionary& dict,
    const polyMesh& mesh
)
:
    tabulatedWallFunction(dict, mesh, typeName),
    interpType_(interpolationTypeNames_[coeffDict_.lookup("interpType")]),
    yPlus_(),
    uPlus_(),
    log10YPlus_(coeffDict_.lookup("log10YPlus")),
    log10UPlus_(coeffDict_.lookup("log10UPlus"))
{
    List<Tuple2<scalar, scalar> > inputTable = coeffDict_.lookup("inputTable");
    if (inputTable.size() < 2)
    {
        FatalErrorIn
        (
            "tabulatedWallFunctions::general::general"
            "("
                "const dictionary&, "
                "const polyMesh&"
            ")"
        )   << "Input table must have at least 2 values" << nl
            << exit(FatalError);
    }

    yPlus_.setSize(inputTable.size());
    uPlus_.setSize(inputTable.size());

    forAll(inputTable, i)
    {
        if (log10YPlus_)
        {
            yPlus_[i] = pow(10, inputTable[i].first());
        }
        else
        {
            yPlus_[i] = inputTable[i].first();
        }

        if (log10UPlus_)
        {
            uPlus_[i] = pow(10, inputTable[i].second());
        }
        else
        {
            uPlus_[i] = inputTable[i].second();
        }
    }

    invertTable();

    if (debug)
    {
        writeData(Info);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tabulatedWallFunctions::general::~general()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::tabulatedWallFunctions::general::yPlus
(
    const scalar uPlus
) const
{
    return interpolate(uPlus, uPlus_, yPlus_);
}


Foam::scalar Foam::tabulatedWallFunctions::general::Re
(
    const scalar uPlus
) const
{
    return uPlus*yPlus(uPlus);
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::tabulatedWallFunctions::general::writeData(Ostream& os) const
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
