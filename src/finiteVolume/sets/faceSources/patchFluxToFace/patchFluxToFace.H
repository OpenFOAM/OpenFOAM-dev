/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::patchFluxToFace

Description
    A topoSetSource to select patch faces according to the flux direction.

Usage
    Example topoSetDict to generate faceSets for inflow and outflow faces
    on the outlet patch:
    \verbatim
    actions
    (
        {
            action  new;
            type    faceSet;
            name    inflow;
            source  patchFluxToFace;
            field   phi;
            patch   outlet;
            inflow  true;
        }

        {
            action  new;
            type    faceSet;
            name    outflow;
            source  patchFluxToFace;
            field   phi;
            patch   outlet;
            inflow  false;
        }
    );
    \endverbatim

SourceFiles
    patchFluxToFace.C

\*---------------------------------------------------------------------------*/

#ifndef patchFluxToFace_H
#define patchFluxToFace_H

#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class patchFluxToFace Declaration
\*---------------------------------------------------------------------------*/

class patchFluxToFace
:
    public topoSetSource
{
    // Private Data

        //- Name of surfaceScalarField
        word fieldName_;

        //- Name of patch
        word patchName_;

        //- Switch for inflow/outflow
        bool inflow_;


    // Private Member Functions

        //- Depending on patchFluxField values add to or delete from faceSet.
        void applyToSet
        (
            const topoSetSource::setAction action,
            const scalarField& patchFluxField,
            topoSet& set
        ) const;


public:

    //- Runtime type information
    TypeName("patchFluxToFace");


    // Constructors

        //- Construct from components
        patchFluxToFace
        (
            const polyMesh& mesh,
            const word& fieldName,
            const word& patchName,
            const bool inflow
        );

        //- Construct from dictionary
        patchFluxToFace
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~patchFluxToFace();


    // Member Functions

        virtual sourceType setType() const
        {
            return FACESETSOURCE;
        }

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet&
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
