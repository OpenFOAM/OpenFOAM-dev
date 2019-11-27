/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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

#include "limitVelocity.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(limitVelocity, 0);
    addToRunTimeSelectionTable
    (
        option,
        limitVelocity,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::limitVelocity::limitVelocity
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    UName_(coeffs_.lookupOrDefault<word>("U", "U")),
    max_(coeffs_.lookup<scalar>("max"))
{
    fieldNames_.setSize(1, UName_);
    applied_.setSize(1, false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::limitVelocity::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.lookup("max") >> max_;

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::fv::limitVelocity::correct(volVectorField& U)
{
    const scalar maxSqrU = sqr(max_);

    vectorField& Uif = U.primitiveFieldRef();

    forAll(cells_, i)
    {
        const label celli = cells_[i];

        const scalar magSqrUi = magSqr(Uif[celli]);

        if (magSqrUi > maxSqrU)
        {
            Uif[celli] *= sqrt(maxSqrU/magSqrUi);
        }
    }

    // handle boundaries in the case of 'all'
    if (selectionMode_ == smAll)
    {
        volVectorField::Boundary& Ubf = U.boundaryFieldRef();

        forAll(Ubf, patchi)
        {
            fvPatchVectorField& Up = Ubf[patchi];

            if (!Up.fixesValue())
            {
                forAll(Up, facei)
                {
                    const scalar magSqrUi = magSqr(Up[facei]);

                    if (magSqrUi > maxSqrU)
                    {
                        Up[facei] *= sqrt(maxSqrU/magSqrUi);
                    }
                }
            }
        }
    }
}


// ************************************************************************* //
