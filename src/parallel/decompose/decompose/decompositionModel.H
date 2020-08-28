/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2020 OpenFOAM Foundation
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
    Foam::decompositionModel

Description
    MeshObject wrapper of decompositionMethod

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef decompositionModel_H
#define decompositionModel_H

#include "IOdictionary.H"
#include "MeshObject.H"
#include "decompositionMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class mapPolyMesh;
class polyMesh;

/*---------------------------------------------------------------------------*\
                     Class decompositionModel Declaration
\*---------------------------------------------------------------------------*/

class decompositionModel
:
    public MeshObject
    <
        polyMesh,
        UpdateableMeshObject,
        decompositionModel
    >,
    public IOdictionary
{

    // Private Data

        mutable autoPtr<decompositionMethod> decomposerPtr_;


public:

    // Declare name of the class and its debug switch
    ClassName("decompositionModel");


    // Selectors

        //- Read (optionally from absolute path) & register on mesh
        static const decompositionModel& New
        (
            const polyMesh& mesh,
            const fileName& decompDictFile = ""
        );

        //- Read (optionally from supplied dictionary) & register on mesh
        static const decompositionModel& New
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const fileName& decompDictFile = ""
        );


    // Constructors

        //- Construct from typeName or optional path to controlDictionary
        decompositionModel(const polyMesh&, const fileName& = "");


        //- Construct from typeName or optional path to controlDictionary
        decompositionModel
        (
            const polyMesh&,
            const dictionary& dict,
            const fileName& = ""
        );


    // Member Functions

        decompositionMethod& decomposer() const
        {
            if (!decomposerPtr_.valid())
            {
                decomposerPtr_ = decompositionMethod::New(*this);
            }
            return decomposerPtr_();
        }

        //- Helper: return IOobject with optionally absolute path provided
        static IOobject selectIO(const IOobject&, const fileName&);


        // UpdateableMeshObject

            virtual bool movePoints()
            {
                return false;
            }

            virtual void updateMesh(const mapPolyMesh&)
            {}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
