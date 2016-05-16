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
    grpUtilitiesFunctionObjects

Description
    This function object evolves a passive scalar transport equation.  The
    field in ininitially zero, to which sources are added.  The field name
    is assigned the name of the function object.  Boundary conditions are
    automatically applied, based on the velocity boundary conditions.

    - the field can be zeroed on start-up using the resetOnStartUp flag
    - to employ the same numerical schemes as the flow velocity, use the
      autoSchemes flag
    - the diffusivity can be set manually using the DT entry, or retrieved
      from the turbulence model (if applicable)

SeeAlso
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    scalarTransport.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_scalarTransport_H
#define functionObjects_scalarTransport_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "surfaceFieldsFwd.H"
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

        //- Name of flux field (optional)
        word phiName_;

        //- Name of velocity field (optional)
        word UName_;

        //- Name of density field (optional)
        word rhoName_;

        //- Diffusion coefficient (optional)
        scalar DT_;

        //- Flag to indicate whether user DT_ is used
        bool userDT_;

        //- Flag to reset scalar field on start-up
        bool resetOnStartUp_;

        //- Number of corrector iterations (optional)
        label nCorr_;

        //- Flag to employ schemes for velocity for scalar transport
        bool autoSchemes_;

        //- Run-time selectable finite volume options, e.g. sources, constraints
        fv::optionList fvOptions_;

        //- The scalar field
        volScalarField T_;


    // Private Member Functions

        //- Return the boundary types for the scalar field
        wordList boundaryTypes() const;

        //- Return the diffusivity field
        tmp<volScalarField> DT(const surfaceScalarField& phi) const;

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
        virtual bool execute(const bool postProcess = false);

        //- Do nothing.
        //  The volScalarField is registered and written automatically
        virtual bool write(const bool postProcess = false);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
