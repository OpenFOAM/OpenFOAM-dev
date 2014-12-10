/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "externalCoupledTemperatureMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "OFstream.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::externalCoupledTemperatureMixedFvPatchScalarField::writeHeader
(
    OFstream& os
) const
{
    os  << "# Values: magSf value qDot htc" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalCoupledTemperatureMixedFvPatchScalarField::
externalCoupledTemperatureMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    externalCoupledMixedFvPatchField<scalar>(p, iF)
{}


Foam::externalCoupledTemperatureMixedFvPatchScalarField::
externalCoupledTemperatureMixedFvPatchScalarField
(
    const externalCoupledTemperatureMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    externalCoupledMixedFvPatchField<scalar>(ptf, p, iF, mapper)
{}


Foam::externalCoupledTemperatureMixedFvPatchScalarField::
externalCoupledTemperatureMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    externalCoupledMixedFvPatchField<scalar>(p, iF, dict)
{}


Foam::externalCoupledTemperatureMixedFvPatchScalarField::
externalCoupledTemperatureMixedFvPatchScalarField
(
    const externalCoupledTemperatureMixedFvPatchScalarField& ecmpf
)
:
    externalCoupledMixedFvPatchField<scalar>(ecmpf)
{}


Foam::externalCoupledTemperatureMixedFvPatchScalarField::
externalCoupledTemperatureMixedFvPatchScalarField
(
    const externalCoupledTemperatureMixedFvPatchScalarField& ecmpf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    externalCoupledMixedFvPatchField<scalar>(ecmpf, iF)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::externalCoupledTemperatureMixedFvPatchScalarField::
~externalCoupledTemperatureMixedFvPatchScalarField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::externalCoupledTemperatureMixedFvPatchScalarField::transferData
(
    OFstream& os
) const
{
    if (log())
    {
        Info<< type() << ": " << this->patch().name()
            << ": writing data to " << os.name()
            << endl;
    }

    const label patchI = patch().index();

    // heat flux [W/m2]
    scalarField qDot(this->patch().size(), 0.0);

    typedef compressible::turbulenceModel cmpTurbModelType;
    static word turbName("turbulenceModel");
    static word thermoName("thermophysicalProperties");

    if (db().foundObject<cmpTurbModelType>(turbName))
    {
        const cmpTurbModelType& turbModel =
            db().lookupObject<cmpTurbModelType>(turbName);

        const basicThermo& thermo = turbModel.thermo();

        const fvPatchScalarField& hep = thermo.he().boundaryField()[patchI];

        qDot = turbModel.alphaEff(patchI)*hep.snGrad();
    }
    else if (db().foundObject<basicThermo>(thermoName))
    {
        const basicThermo& thermo = db().lookupObject<basicThermo>(thermoName);

        const fvPatchScalarField& hep = thermo.he().boundaryField()[patchI];

        qDot = thermo.alpha().boundaryField()[patchI]*hep.snGrad();
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::externalCoupledTemperatureMixedFvPatchScalarField::"
            "transferData"
            "("
                "OFstream&"
            ") const"
        )   << "Condition requires either compressible turbulence and/or "
            << "thermo model to be available" << exit(FatalError);
    }

    // patch temperature [K]
    const scalarField Tp(*this);

    // near wall cell temperature [K]
    const scalarField Tc(patchInternalField());

    // heat transfer coefficient [W/m2/K]
    const scalarField htc(qDot/(Tp - Tc + ROOTVSMALL));

    if (Pstream::parRun())
    {
        int tag = Pstream::msgType() + 1;

        List<Field<scalar> > magSfs(Pstream::nProcs());
        magSfs[Pstream::myProcNo()].setSize(this->patch().size());
        magSfs[Pstream::myProcNo()] = this->patch().magSf();
        Pstream::gatherList(magSfs, tag);

        List<Field<scalar> > values(Pstream::nProcs());
        values[Pstream::myProcNo()].setSize(this->patch().size());
        values[Pstream::myProcNo()] = Tp;
        Pstream::gatherList(values, tag);

        List<Field<scalar> > qDots(Pstream::nProcs());
        qDots[Pstream::myProcNo()].setSize(this->patch().size());
        qDots[Pstream::myProcNo()] = qDot;
        Pstream::gatherList(qDots, tag);

        List<Field<scalar> > htcs(Pstream::nProcs());
        htcs[Pstream::myProcNo()].setSize(this->patch().size());
        htcs[Pstream::myProcNo()] = htc;
        Pstream::gatherList(htcs, tag);

        if (Pstream::master())
        {
            forAll(values, procI)
            {
                const Field<scalar>& magSf = magSfs[procI];
                const Field<scalar>& value = values[procI];
                const Field<scalar>& qDot = qDots[procI];
                const Field<scalar>& htc = htcs[procI];

                forAll(magSf, faceI)
                {
                    os  << magSf[faceI] << token::SPACE
                        << value[faceI] << token::SPACE
                        << qDot[faceI] << token::SPACE
                        << htc[faceI] << token::SPACE
                        << nl;
                }
            }

            os.flush();
        }
    }
    else
    {
        const Field<scalar>& magSf(this->patch().magSf());

        forAll(patch(), faceI)
        {
            os  << magSf[faceI] << token::SPACE
                << Tp[faceI] << token::SPACE
                << qDot[faceI] << token::SPACE
                << htc[faceI] << token::SPACE
                << nl;
        }

        os.flush();
    }
}


void Foam::externalCoupledTemperatureMixedFvPatchScalarField::evaluate
(
    const Pstream::commsTypes comms
)
{
    externalCoupledMixedFvPatchField<scalar>::evaluate(comms);
}


void Foam::externalCoupledTemperatureMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    externalCoupledMixedFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        externalCoupledTemperatureMixedFvPatchScalarField
    );
}


// ************************************************************************* //
