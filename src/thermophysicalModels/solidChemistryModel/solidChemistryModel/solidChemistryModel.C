/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "solidChemistryModel.H"
#include "reactingMixture.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class SolidThermo>
Foam::solidChemistryModel<CompType, SolidThermo>::
solidChemistryModel
(
    const fvMesh& mesh
)
:
    CompType(mesh),
    ODESystem(),
    Ys_(this->solidThermo().composition().Y()),
    reactions_
    (
        dynamic_cast<const reactingMixture<SolidThermo>& >
        (
            this->solidThermo()
        )
    ),
    solidThermo_
    (
        dynamic_cast<const reactingMixture<SolidThermo>& >
        (
            this->solidThermo()
        ).speciesData()
    ),
    nSolids_(Ys_.size()),
    nReaction_(reactions_.size()),
    RRs_(nSolids_),
    reactingCells_(mesh.nCells(), true)
{
    // create the fields for the chemistry sources
    forAll(RRs_, fieldI)
    {
        RRs_.set
        (
            fieldI,
            new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    "RRs." + Ys_[fieldI].name(),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
            )
        );
   }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class SolidThermo>
Foam::solidChemistryModel<CompType, SolidThermo>::
~solidChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class SolidThermo>
Foam::scalar Foam::solidChemistryModel<CompType, SolidThermo>::solve
(
    const scalarField& deltaT
)
{
    notImplemented
    (
        "solidChemistryModel::solve(const scalarField& deltaT)"
    );
    return 0;
}


template<class CompType, class SolidThermo>
Foam::tmp<Foam::volScalarField>
Foam::solidChemistryModel<CompType, SolidThermo>::tc() const
{
    notImplemented
    (
        "solidChemistryModel::tc()"
    );
    return volScalarField::null();
}


template<class CompType, class SolidThermo>
Foam::tmp<Foam::volScalarField>
Foam::solidChemistryModel<CompType, SolidThermo>::Sh() const
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                "Sh",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        scalarField& Sh = tSh();

        forAll(Ys_, i)
        {
            forAll(Sh, cellI)
            {
                scalar hf = solidThermo_[i].Hc();
                Sh[cellI] -= hf*RRs_[i][cellI];
            }
        }
    }

    return tSh;
}


template<class CompType, class SolidThermo>
Foam::tmp<Foam::volScalarField>
Foam::solidChemistryModel<CompType, SolidThermo>::dQ() const
{
    tmp<volScalarField> tdQ
    (
        new volScalarField
        (
            IOobject
            (
                "dQ",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("dQ", dimEnergy/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        volScalarField& dQ = tdQ();
        dQ.dimensionedInternalField() = this->mesh_.V()*Sh()();
    }

    return tdQ;
}


template<class CompType, class SolidThermo>
void Foam::solidChemistryModel<CompType, SolidThermo>::setCellReacting
(
    const label cellI,
    const bool active
)
{
    reactingCells_[cellI] = active;
}

// ************************************************************************* //
