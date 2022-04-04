/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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
    Foam::fvMeshMovers::inkJet

Description
    Mesh motion specifically for the "pumping" system of an ink-jet
    injector.

    The set of points in the "pumping" region are compressed and expanded
    sinusoidally to impose a sinusoidal variation of the flow at the
    nozzle exit.

SourceFiles
    inkJetFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshMoversInkJet_H
#define fvMeshMoversInkJet_H

#include "fvMeshMover.H"
#include "dictionary.H"
#include "pointIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshMovers
{

/*---------------------------------------------------------------------------*\
                         Class inkJet Declaration
\*---------------------------------------------------------------------------*/

class inkJet
:
    public fvMeshMover
{
    // Private Data

        dictionary meshCoeffs_;

        scalar amplitude_;
        scalar frequency_;
        scalar refPlaneX_;

        pointIOField stationaryPoints_;

        //- Optional list of vectorFields to update for mesh motion
        //  For modern solvers using Uf and correctPhi to update the flux
        //  after motion it is not necessary to specify a "velocityFields" list
        velocityMotionCorrection velocityMotionCorrection_;


public:

    //- Runtime type information
    TypeName("inkJet");


    // Constructors

        //- Construct from fvMesh
        inkJet(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        inkJet(const inkJet&) = delete;


    //- Destructor
    ~inkJet();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const inkJet&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshMovers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
