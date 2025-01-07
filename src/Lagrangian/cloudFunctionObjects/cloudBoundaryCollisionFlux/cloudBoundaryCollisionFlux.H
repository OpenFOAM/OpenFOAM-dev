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

Class
    Foam::functionObjects::cloudBoundaryCollisionFlux

Description
    Base class for functions which generate a boundary collision flux for a
    cloud

SourceFiles
    cloudBoundaryCollisionFlux.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_cloudBoundaryCollisionFlux_H
#define functionObjects_cloudBoundaryCollisionFlux_H

#include "LagrangianPatchField.H"
#include "cloudFvMeshFunctionObject.H"
#include "surfaceFields.H"
#include "LagrangianFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                 Class cloudBoundaryCollisionFlux Declaration
\*---------------------------------------------------------------------------*/

class cloudBoundaryCollisionFlux
:
    public cloudFvMeshFunctionObject
{
private:

    // Private Data

        //- Field of face indices
        autoPtr<LagrangianLabelField> faceiPtr_;

        //- Field of the property for which to compute the flux
        autoPtr<LagrangianScalarField> qPtr_;

        //- Name of the flux
        const word phiName_;

        //- Dimensions of the flux
        const dimensionSet phiDims_;

        //- Flux boundary field
        surfaceScalarField::Boundary phib_;


    // Private Member Functions

        //- Return the property for which to compute the flux
        virtual tmp<LagrangianSubScalarSubField> q
        (
            const LagrangianSubScalarSubField& fraction,
            const label sign
        ) const = 0;


public:

    //- Runtime type information
    TypeName("cloudBoundaryCollisionFlux");


    // Constructors

        //- Construct from Time and dictionary
        cloudBoundaryCollisionFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const word& phiName,
            const dimensionSet& phiDims
        );

        //- Disallow default bitwise copy construction
        cloudBoundaryCollisionFlux(const cloudBoundaryCollisionFlux&) = delete;


    //- Destructor
    virtual ~cloudBoundaryCollisionFlux();


    // Member Functions

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Return false so this function does not executed at the start
        virtual bool executeAtStart() const;

        //- Calculate and store the number flux
        virtual bool execute();

        //- Hook before solution steps
        virtual void preSolve();

        //- Hook before all face crossings
        virtual void preCrossFaces
        (
            const LagrangianScalarInternalDynamicField& fraction
        );

        //- Hook before face crossings of a specific sub-mesh
        virtual void preCrossFaces
        (
            const LagrangianSubScalarSubField& fraction
        );

        //- Hook following face crossings of a specific sub-mesh
        virtual void postCrossFaces
        (
            const LagrangianSubScalarSubField& fraction
        );

        //- Hook after all face crossings
        virtual void postCrossFaces
        (
            const LagrangianScalarInternalDynamicField& fraction
        );

        //- Write the number flux
        virtual bool write();

        //- Clear the number flux
        virtual bool clear();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cloudBoundaryCollisionFlux&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
