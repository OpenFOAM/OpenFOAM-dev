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
    Foam::LagrangianInjection

Description
    Base class for Lagrangian injections. Minimal wrapper over
    LagrangianSource. Implements some utility functions, and reports that all
    fields require a source value to be specified. Also serves as a type that
    can be identified and cast to.

SourceFiles
    LagrangianInjection.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianInjection_H
#define LagrangianInjection_H

#include "LagrangianModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class timeIOdictionary;

/*---------------------------------------------------------------------------*\
                    Class LagrangianInjection Declaration
\*---------------------------------------------------------------------------*/

class LagrangianInjection
:
    public LagrangianModel
{
protected:

    // Protected Member Functions

        //- Check the location and error or warn as appropriate
        void checkLocation
        (
            const LagrangianMesh::location,
            const point&
        ) const;

        //- Check the locations and error or warn as appropriate
        void checkLocation
        (
            const List<LagrangianMesh::location>&,
            const List<point>&
        ) const;

        //- Filter out injected elements that are on other processes. These
        //  are identified by cell/face/faceTri indices equal to -1.
        template<class ... OtherFields>
        static void filter
        (
            barycentricField& coordinates,
            labelField& celli,
            labelField& facei,
            labelField& faceTrii,
            OtherFields& ... otherFields
        );


public:

    //- Runtime type information
    TypeName("LagrangianInjection");


    // Constructors

        //- Construct from components
        LagrangianInjection(const word& name, const LagrangianMesh& mesh);


    //- Destructor
    virtual ~LagrangianInjection();


    // Member Functions

        // Checks

            //- Return true if the LagrangianModel adds a source term to the
            //  given field's transport equation. Returns false. Injections are
            //  instantaneous. They do not create continuous sources.
            virtual bool addsSupToField(const word& fieldName) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LagrangianInjectionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
