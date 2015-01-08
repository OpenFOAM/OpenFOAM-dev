/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::wallDist

Description
    Interface to run-time selectable methods to calculate the distance-to-wall
    field.

SourceFiles
    wallDist.C

\*---------------------------------------------------------------------------*/

#ifndef wallDist_H
#define wallDist_H

#include "MeshObject.H"
#include "patchDistMethod.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class wallDist Declaration
\*---------------------------------------------------------------------------*/

class wallDist
:
    public MeshObject<fvMesh, UpdateableMeshObject, wallDist>
{
    // Private data

        //- Run-time selected method to generate the distance-to-wall field
        autoPtr<patchDistMethod> pdm_;

        //- Distance-to-wall field
        volScalarField y_;

        //- Distance-to-wall field
        tmp<volVectorField> n_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        wallDist(const wallDist&);

        //- Disallow default bitwise assignment
        void operator=(const wallDist&);


public:

    // Declare name of the class and its debug switch
    ClassName("wallDist");


    // Constructors

        //- Construct from mesh
        wallDist(const fvMesh& mesh);


    //- Destructor
    virtual ~wallDist();


    // Member Functions

        //- Return reference to cached distance-to-wall field
        const volScalarField& y() const
        {
            return y_;
        }

        //- Return reference to cached normal-to-wall field
        const volVectorField& n() const
        {
            return n_();
        }

        //- Update the y-field when the mesh moves
        virtual bool movePoints();

        //- Update the y-field when the mesh changes
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
