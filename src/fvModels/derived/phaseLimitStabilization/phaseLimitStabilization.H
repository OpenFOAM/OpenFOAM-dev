/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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
    Foam::fv::phaseLimitStabilization

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
        type            phaseLimitStabilization;

        field           sigma.liquid;

        rate            rLambda.liquid;
        residualAlpha   1e-3;
    }
    \endverbatim

SourceFiles
    phaseLimitStabilization.C

\*---------------------------------------------------------------------------*/

#ifndef phaseLimitStabilization_H
#define phaseLimitStabilization_H

#include "fvModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                   Class phaseLimitStabilization Declaration
\*---------------------------------------------------------------------------*/

class phaseLimitStabilization
:
    public fvModel
{
    // Private Data

        //- Field name
        word fieldName_;

        //- Rate field name
        word rateName_;

        //- Residual alpha value below which stabilization is applied
        scalar residualAlpha_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();

        //- Add source terms to a phase equation
        template<class Type>
        void addSupType
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            fvMatrix<Type>& eqn,
            const word& fieldName
        ) const;


public:

    //- Runtime type information
    TypeName("phaseLimitStabilization");


    // Constructors

        //- Construct from components
        phaseLimitStabilization
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        phaseLimitStabilization(const phaseLimitStabilization&) = delete;


    //- Destructor
    virtual ~phaseLimitStabilization()
    {}


    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Sources

            //- Add a source term to a phase equation
            FOR_ALL_FIELD_TYPES(DEFINE_FV_MODEL_ADD_ALPHA_RHO_SUP);


        //- Read dictionary
        virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const phaseLimitStabilization&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
