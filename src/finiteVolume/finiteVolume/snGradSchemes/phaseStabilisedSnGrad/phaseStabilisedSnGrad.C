/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "phaseStabilisedSnGrad.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "localMin.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
phaseStabilisedSnGrad<Type>::phaseStabilisedSnGrad
(
    const fvMesh& mesh,
    Istream& schemeData
)
:
    snGradScheme<Type>(mesh),
    correctedScheme_
    (
        fv::snGradScheme<Type>::New(this->mesh(), schemeData)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
phaseStabilisedSnGrad<Type>::~phaseStabilisedSnGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<surfaceScalarField> phaseStabilisedSnGrad<Type>::deltaCoeffs
(
    const VolField<Type>& vf
) const
{
    return correctedScheme_->deltaCoeffs(vf);
}


template<class Type>
tmp<SurfaceField<Type>>
phaseStabilisedSnGrad<Type>::correction
(
    const VolField<Type>& vf
) const
{
    const SurfaceField<Type> corr
    (
        correctedScheme_().correction(vf)
    );

    const volScalarField& alpha
    (
        vf.db().template lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", vf.group())
        )
    );

    const surfaceScalarField limiter
    (
        pos(localMin<scalar>(vf.mesh()).interpolate(alpha) - 1e-3)
    );

    if (fv::debug)
    {
        InfoInFunction
            << "limiter min: " << min(limiter.primitiveField())
            << " max: "<< max(limiter.primitiveField())
            << " avg: " << average(limiter.primitiveField()) << endl;
    }

    return limiter*corr;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
