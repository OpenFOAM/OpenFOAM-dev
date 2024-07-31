/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2024 OpenFOAM Foundation
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

#include "limitMag.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(limitMag, 0);
    addToRunTimeSelectionTable
    (
        fvConstraint,
        limitMag,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::limitMag::readCoeffs(const dictionary& dict)
{
    fieldName_ = dict.lookup<word>("field");
    max_ = dict.lookup<scalar>("max");
}


template<class Type>
inline bool Foam::fv::limitMag::constrainType
(
    VolField<Type>& psi
) const
{
    const scalar maxSqrPsi = sqr(max_);

    Field<Type>& psiif = psi.primitiveFieldRef();

    const labelUList cells = set_.cells();

    forAll(cells, i)
    {
        const label celli = cells[i];

        const scalar magSqrPsii = magSqr(psiif[celli]);

        if (magSqrPsii > maxSqrPsi)
        {
            psiif[celli] *= sqrt(maxSqrPsi/magSqrPsii);
        }
    }

    // handle boundaries in the case of 'all'
    if (set_.selectionType() == fvCellSet::selectionTypes::all)
    {
        typename VolField<Type>::Boundary& psibf =
            psi.boundaryFieldRef();

        forAll(psibf, patchi)
        {
            fvPatchField<Type>& psip = psibf[patchi];

            if (!psip.fixesValue())
            {
                forAll(psip, facei)
                {
                    const scalar magSqrPsii = magSqr(psip[facei]);

                    if (magSqrPsii > maxSqrPsi)
                    {
                        psip[facei] *= sqrt(maxSqrPsi/magSqrPsii);
                    }
                }
            }
        }
    }

    return cells.size();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::limitMag::limitMag
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvConstraint(name, modelType, mesh, dict),
    set_(mesh, coeffs(dict)),
    fieldName_(word::null),
    max_(vGreat)
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::limitMag::constrainedFields() const
{
    return wordList(1, fieldName_);
}


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_CONSTRAINT_CONSTRAIN_FIELD,
    fv::limitMag
);


bool Foam::fv::limitMag::movePoints()
{
    set_.movePoints();
    return true;
}


void Foam::fv::limitMag::topoChange(const polyTopoChangeMap& map)
{
    set_.topoChange(map);
}


void Foam::fv::limitMag::mapMesh(const polyMeshMap& map)
{
    set_.mapMesh(map);
}


void Foam::fv::limitMag::distribute(const polyDistributionMap& map)
{
    set_.distribute(map);
}


bool Foam::fv::limitMag::read(const dictionary& dict)
{
    if (fvConstraint::read(dict))
    {
        set_.read(coeffs(dict));
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
