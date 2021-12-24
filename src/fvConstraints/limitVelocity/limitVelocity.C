/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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
        fvConstraint,
        limitVelocity,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::limitVelocity::readCoeffs()
{
    UName_ = coeffs().lookupOrDefault<word>("U", "U");
    max_ = coeffs().lookup<scalar>("max");
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
    fvConstraint(name, modelType, dict, mesh),
    set_(coeffs(), mesh),
    UName_(word::null),
    max_(vGreat)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::limitVelocity::constrainedFields() const
{
    return wordList(1, UName_);
}


bool Foam::fv::limitVelocity::constrain(volVectorField& U) const
{
    const scalar maxSqrU = sqr(max_);

    vectorField& Uif = U.primitiveFieldRef();

    const labelList& cells = set_.cells();

    forAll(cells, i)
    {
        const label celli = cells[i];

        const scalar magSqrUi = magSqr(Uif[celli]);

        if (magSqrUi > maxSqrU)
        {
            Uif[celli] *= sqrt(maxSqrU/magSqrUi);
        }
    }

    // handle boundaries in the case of 'all'
    if (set_.selectionMode() == fvCellSet::selectionModeType::all)
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

    return cells.size();
}


void Foam::fv::limitVelocity::updateMesh(const mapPolyMesh& map)
{
    set_.updateMesh(map);
}


void Foam::fv::limitVelocity::distribute(const mapDistributePolyMesh& map)
{
    set_.distribute(map);
}


bool Foam::fv::limitVelocity::movePoints()
{
    set_.movePoints();
    return true;
}


bool Foam::fv::limitVelocity::read(const dictionary& dict)
{
    if (fvConstraint::read(dict))
    {
        set_.read(coeffs());
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
