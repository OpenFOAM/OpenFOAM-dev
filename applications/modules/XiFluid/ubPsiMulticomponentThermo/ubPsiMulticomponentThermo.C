/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "ubPsiMulticomponentThermo.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ubPsiMulticomponentThermo, 0);
    defineRunTimeSelectionTable(ubPsiMulticomponentThermo, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ubPsiMulticomponentThermo::implementation::implementation
(
    const dictionary& dict,
    const wordList& specieNames,
    const fvMesh& mesh,
    const word& phaseName
)
:
    species_(specieNames),
    Y_(species_.size())
{
    forAll(species_, i)
    {
        Y_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName(species_[i], phaseName),
                    mesh.time().name(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::ubPsiMulticomponentThermo>
Foam::ubPsiMulticomponentThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<ubPsiMulticomponentThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ubPsiMulticomponentThermo::~ubPsiMulticomponentThermo()
{}


Foam::ubPsiMulticomponentThermo::implementation::~implementation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::speciesTable&
Foam::ubPsiMulticomponentThermo::implementation::species() const
{
    return species_;
}


Foam::PtrList<Foam::volScalarField>&
Foam::ubPsiMulticomponentThermo::implementation::Y()
{
    return Y_;
}


const Foam::PtrList<Foam::volScalarField>&
Foam::ubPsiMulticomponentThermo::implementation::Y() const
{
    return Y_;
}


// ************************************************************************* //
