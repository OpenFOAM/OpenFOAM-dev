/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2023 OpenFOAM Foundation
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
    Foam::singleRegionSolutionControl

Description
    Single-region-specific derivation of the solution control class

SourceFiles
    singleRegionSolutionControlI.H
    singleRegionSolutionControl.C
    singleRegionSolutionControlTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef singleRegionSolutionControl_H
#define singleRegionSolutionControl_H

#include "fvMesh.H"
#include "solutionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class singleRegionSolutionControl Declaration
\*---------------------------------------------------------------------------*/

class singleRegionSolutionControl
:
    public solutionControl
{
protected:

    // Protected data

        //- Reference to the mesh
        fvMesh& mesh_;


    // Protected Member Functions

        //- Return true if the dependencies have been modified
        virtual bool dependenciesModified() const;

        //- Store previous iteration for fields
        template<class Type>
        void storePrevIterTypeFields() const;


public:

    // Static Data Members

        //- Run-time type information
        TypeName("singleRegionSolutionControl");


    // Constructors

        //- Construct from a mesh and the name of the algorithm
        singleRegionSolutionControl(fvMesh& mesh, const word& algorithmName);


    //- Destructor
    virtual ~singleRegionSolutionControl();


    // Member Functions

        // Access

            //- Return the mesh
            inline const fvMesh& mesh() const;

            //- Return the mesh
            inline fvMesh& mesh();

            //- Return the solution dictionary
            virtual const dictionary& dict() const;


        // Evolution

            //- Update the mesh data for the final state
            void updateFinal(const bool finalIter) const;

            //- Store previous iteration for fields
            void storePrevIterFields();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "singleRegionSolutionControlTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "singleRegionSolutionControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
