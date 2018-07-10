/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

Class
    Foam::fv::PhaseLimitStabilization

Description
    Stabilization source for phase transport equations

    Applies an implicit source to the phase transport equation for the specified
    \c field when the phase volume fraction is below \c residualAlpha.  The
    stabilization rate is provided by the registered
    uniformDimensionedScalarField \c rate, which could be extended to also
    support volScalarField and volScalarField::Internal field types.  The \c
    field is currently stabilized towards zero in the limit of the phase volume
    fraction approaching zero but this could be extended to support a
    specified value or a value or field looked-up from the database.

Usage
    Example usage:
    \verbatim
    stabilization
    {
        type            symmTensorPhaseLimitStabilization;

        field           sigma.liquid;
        rate            rLambda.liquid;
        residualAlpha   1e-3;
    }
    \endverbatim

SourceFiles
    PhaseLimitStabilization.C

\*---------------------------------------------------------------------------*/

#ifndef PhaseLimitStabilization_H
#define PhaseLimitStabilization_H

#include "fvOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                Class PhaseLimitStabilization Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class PhaseLimitStabilization
:
    public option
{
    // Private data

        //- Field name
        word fieldName_;

        //- Rate field name
        word rateName_;

        //- Residual alpha value below which stabilization is applied
        scalar residualAlpha_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        PhaseLimitStabilization(const PhaseLimitStabilization&);

        //- Disallow default bitwise assignment
        void operator=(const PhaseLimitStabilization&);


public:

    //- Runtime type information
    TypeName("PhaseLimitStabilization");


    // Constructors

        //- Construct from components
        PhaseLimitStabilization
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~PhaseLimitStabilization()
    {}


    // Member Functions

        //- Source term to compressible phase equation
        virtual void addSup
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            fvMatrix<Type>& eqn,
            const label fieldi
        );

        //- Read dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PhaseLimitStabilization.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
