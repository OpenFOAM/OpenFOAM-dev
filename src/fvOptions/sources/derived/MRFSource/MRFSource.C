/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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

#include "MRFSource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "MRFZone.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(MRFSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        MRFSource,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::fv::MRFSource::initialise()
{
    if (selectionMode_ != smCellZone)
    {
        FatalErrorIn("void Foam::MRFSource::initialise()")
            << "The MRF region must be specified as a cellZone. Current "
            << "selection mode is " << selectionModeTypeNames_[selectionMode_]
            << exit(FatalError);
    }

    mrfPtr_.reset
    (
        new MRFZone
        (
            name_,
            mesh_,
            coeffs_,
            cellSetName_
        )
    );

    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);

    mrfPtr_->correctBoundaryVelocity(const_cast<volVectorField&>(U));

    fieldNames_.setSize(1, UName_);
    applied_.setSize(1, false);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::MRFSource::MRFSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    mrfPtr_(NULL),
    UName_(coeffs_.lookupOrDefault<word>("UName", "U"))
{
    initialise();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::MRFSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    // Add to rhs of equation
    mrfPtr_->addCoriolis(eqn, true);
}


void Foam::fv::MRFSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    // Add to rhs of equation
    mrfPtr_->addCoriolis(rho, eqn, true);
}


void Foam::fv::MRFSource::makeRelative(surfaceScalarField& phi) const
{
    mrfPtr_->makeRelative(phi);
}


void Foam::fv::MRFSource::makeRelative
(
    FieldField<fvsPatchField, scalar>& phi
) const
{
    mrfPtr_->makeRelative(phi);
}


void Foam::fv::MRFSource::makeRelative
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    mrfPtr_->makeRelative(rho, phi);
}


void Foam::fv::MRFSource::makeAbsolute(surfaceScalarField& phi) const
{
    mrfPtr_->makeAbsolute(phi);
}


void Foam::fv::MRFSource::makeAbsolute
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    mrfPtr_->makeAbsolute(rho, phi);
}


void Foam::fv::MRFSource::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::fv::MRFSource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        coeffs_.readIfPresent("UName", UName_);

        initialise();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
