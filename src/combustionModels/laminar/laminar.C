/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "laminar.H"
#include "fvmSup.H"
#include "localEulerDdtScheme.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::combustionModels::laminar<Type>::laminar
(
    const word& modelType,
    const fvMesh& mesh
)
:
    Type(modelType, mesh),
    integrateReactionRate_
    (
        this->coeffs().lookupOrDefault("integrateReactionRate", true)
    )
{
    if (integrateReactionRate_)
    {
        Info<< "    using integrated reaction rate" << endl;
    }
    else
    {
        Info<< "    using instantaneous reaction rate" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::combustionModels::laminar<Type>::~laminar()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::laminar<Type>::tc() const
{
    return this->chemistryPtr_->tc();
}


template<class Type>
void Foam::combustionModels::laminar<Type>::correct()
{
    if (this->active())
    {
        if (integrateReactionRate_)
        {
            word ddtScheme(this->mesh().ddtScheme("Yi"));

            if (ddtScheme == fv::localEulerDdtScheme<scalar>::typeName)
            {
                const scalarField& rDeltaT =
                    this->mesh().objectRegistry::
                    template lookupObject<volScalarField>
                    (
                        "rDeltaT"
                    );

                if (this->coeffs().found("maxIntegrationTime"))
                {
                    scalar maxIntegrationTime
                    (
                        readScalar(this->coeffs().lookup("maxIntegrationTime"))
                    );

                    this->chemistryPtr_->solve
                    (
                        min(1.0/rDeltaT, maxIntegrationTime)()
                    );
                }
                else
                {
                    this->chemistryPtr_->solve((1.0/rDeltaT)());
                }
            }
            else
            {
                this->chemistryPtr_->solve(this->mesh().time().deltaTValue());
            }
        }
        else
        {
            this->chemistryPtr_->calculate();
        }
    }
}


template<class Type>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::laminar<Type>::R(volScalarField& Y) const
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Y, dimMass/dimTime));

    fvScalarMatrix& Su = tSu();

    if (this->active())
    {
        const label specieI = this->thermo().composition().species()[Y.name()];

        Su += this->chemistryPtr_->RR(specieI);
    }

    return tSu;
}


template<class Type>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::laminar<Type>::dQ() const
{
    tmp<volScalarField> tdQ
    (
        new volScalarField
        (
            IOobject
            (
                typeName + ":dQ",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("dQ", dimEnergy/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->active())
    {
        tdQ() = this->chemistryPtr_->dQ();
    }

    return tdQ;
}


template<class Type>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::laminar<Type>::Sh() const
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                typeName + ":Sh",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->active())
    {
        tSh() = this->chemistryPtr_->Sh();
    }

    return tSh;
}


template<class Type>
bool Foam::combustionModels::laminar<Type>::read()
{
    if (Type::read())
    {
        this->coeffs().lookup("integrateReactionRate")
            >> integrateReactionRate_;
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
