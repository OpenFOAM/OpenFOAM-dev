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
    Foam::cloudFunctionObjectUList

Description
    List of references to the cloud function objects. Designed to be
    constructed temporarily for the scope of cloud evolution so that
    cloud-specific hooks can be used.

SourceFiles
    cloudFunctionObjectUList.C

\*---------------------------------------------------------------------------*/

#ifndef cloudFunctionObjectUList_H
#define cloudFunctionObjectUList_H

#include "cloudFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class cloudFunctionObjectUList Declaration
\*---------------------------------------------------------------------------*/

class cloudFunctionObjectUList
:
    public UPtrList<functionObjects::cloudFunctionObject>
{
private:

    // Private Data

        //- Is this in an inner scope? Prevents calls to pre and post solve
        const bool inner_;


public:

    //- Runtime type information
    TypeName("cloudFunctionObjectUList");


    // Constructors

        //- Construct from a cloud
        cloudFunctionObjectUList(const cloud&, const bool inner = false);

        //- Disallow default bitwise copy construction
        cloudFunctionObjectUList
        (
            const cloudFunctionObjectUList&
        ) = delete;


    //- Destructor
    virtual ~cloudFunctionObjectUList();


    // Member Functions

        //- Solve equations and/or update continually changing properties
        virtual void calculate
        (
            const LagrangianSubScalarField& deltaT,
            const bool final
        );

        //- Hook before face crossings
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

        //- Hook after face crossings
        virtual void postCrossFaces
        (
            const LagrangianScalarInternalDynamicField& fraction
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cloudFunctionObjectUList&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
