/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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
    Foam::functionObjects::scalarTransport

Group
    grpSolversFunctionObjects

Description
    This function object evolves a passive scalar transport equation.

    - To specify the field name set the 'field' entry
    - To employ the same numerical schemes as another field set
      the 'schemesField' entry,
    - The diffusivity can be set manually using the 'D' entry, or retrieved
      from the turbulence model (if applicable).

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    scalarTransport.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_scalarTransport_H
#define functionObjects_scalarTransport_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "fvOptionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class scalarTransport Declaration
\*---------------------------------------------------------------------------*/

class scalarTransport
:
    public fvMeshFunctionObject
{
    // Private data

        //- Name of field to process
        word fieldName_;

        //- Name of flux field (optional)
        word phiName_;

        //- Name of density field (optional)
        word rhoName_;

        //- Diffusion coefficient (optional)
        scalar D_;

        //- Flag to indicate whether a constant, uniform D_ is specified
        bool constantD_;

        //- Number of corrector iterations (optional)
        label nCorr_;

        //- Name of field whose schemes are used (optional)
        word schemesField_;

        //- Run-time selectable finite volume options, e.g. sources, constraints
        fv::optionList fvOptions_;

        //- The scalar field
        volScalarField s_;


    // Private Member Functions

        //- Return the diffusivity field
        tmp<volScalarField> D(const surfaceScalarField& phi) const;

        //- Disallow default bitwise copy construct
        scalarTransport(const scalarTransport&);

        //- Disallow default bitwise assignment
        void operator=(const scalarTransport&);


public:

    //- Runtime type information
    TypeName("scalarTransport");


    // Constructors

        //- Construct from Time and dictionary
        scalarTransport
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~scalarTransport();


    // Member Functions

        //- Read the scalarTransport data
        virtual bool read(const dictionary&);

        //- Calculate the scalarTransport
        virtual bool execute();

        //- Do nothing.
        //  The volScalarField is registered and written automatically
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
