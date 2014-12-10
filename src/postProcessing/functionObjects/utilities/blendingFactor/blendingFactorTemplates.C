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

#include "gaussConvectionScheme.H"
#include "blendedSchemeBase.H"
#include "fvcCellReduce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::volScalarField& Foam::blendingFactor::factor
(
    const GeometricField<Type, fvPatchField, volMesh>& field
)
{
    const word fieldName = "blendingFactor:" + field.name();

    if (!obr_.foundObject<volScalarField>(fieldName))
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        volScalarField* factorPtr =
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("0", dimless, 0.0),
                zeroGradientFvPatchScalarField::typeName
            );

        obr_.store(factorPtr);
    }

    return
        const_cast<volScalarField&>
        (
            obr_.lookupObject<volScalarField>(fieldName)
        );
}


template<class Type>
void Foam::blendingFactor::calc()
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (!obr_.foundObject<fieldType>(fieldName_))
    {
        return;
    }

    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    const fieldType& field = mesh.lookupObject<fieldType>(fieldName_);

    const word divScheme("div(" + phiName_ + ',' + fieldName_ + ')');
    ITstream& its = mesh.divScheme(divScheme);

    const surfaceScalarField& phi =
        mesh.lookupObject<surfaceScalarField>(phiName_);

    tmp<fv::convectionScheme<Type> > cs =
        fv::convectionScheme<Type>::New(mesh, phi, its);

    const fv::gaussConvectionScheme<Type>& gcs =
        refCast<const fv::gaussConvectionScheme<Type> >(cs());

    const surfaceInterpolationScheme<Type>& interpScheme =
        gcs.interpScheme();

    if (!isA<blendedSchemeBase<Type> >(interpScheme))
    {
        FatalErrorIn("void Foam::blendingFactor::execute()")
            << interpScheme.typeName << " is not a blended scheme"
            << exit(FatalError);
    }

    // retrieve the face-based blending factor
    const blendedSchemeBase<Type>& blendedScheme =
        refCast<const blendedSchemeBase<Type> >(interpScheme);
    const surfaceScalarField factorf(blendedScheme.blendingFactor(field));

    // convert into vol field whose values represent the local face maxima
    volScalarField& factor = this->factor(field);
    factor = fvc::cellReduce(factorf, maxEqOp<scalar>());
    factor.correctBoundaryConditions();
}


// ************************************************************************* //
