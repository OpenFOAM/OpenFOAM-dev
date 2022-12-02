/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "limitedSnGrad.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
limitedSnGrad<Type>::limitedSnGrad(const fvMesh& mesh)
:
    snGradScheme<Type>(mesh),
    correctedScheme_(new correctedSnGrad<Type>(this->mesh())),
    limitCoeff_(1)
{}


template<class Type>
limitedSnGrad<Type>::limitedSnGrad(const fvMesh& mesh, Istream& schemeData)
:
    snGradScheme<Type>(mesh),
    correctedScheme_(lookupCorrectedScheme(schemeData))
{
    if (limitCoeff_ < 0 || limitCoeff_ > 1)
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "limitCoeff is specified as " << limitCoeff_
            << " but should be >= 0 && <= 1"
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
limitedSnGrad<Type>::~limitedSnGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<surfaceScalarField> limitedSnGrad<Type>::deltaCoeffs
(
    const VolField<Type>& vf
) const
{
    return correctedScheme_->deltaCoeffs(vf);
}


template<class Type>
tmp<SurfaceField<Type>>
limitedSnGrad<Type>::correction
(
    const VolField<Type>& vf
) const
{
    const SurfaceField<Type> corr
    (
        correctedScheme_().correction(vf)
    );

    const surfaceScalarField limiter
    (
        min
        (
            limitCoeff_
           *mag(snGradScheme<Type>::snGrad(vf, deltaCoeffs(vf), "SndGrad"))
           /(
                (1 - limitCoeff_)*mag(corr)
              + dimensionedScalar(corr.dimensions(), small)
            ),
            dimensionedScalar(dimless, 1.0)
        )
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
