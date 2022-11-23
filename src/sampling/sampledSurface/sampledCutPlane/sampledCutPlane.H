/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::sampledSurfaces::cutPlane

Description
    A sampledSurface defined by a plane.

    Example:
    \verbatim
    {
        type        cutPlane;
        planeType   pointAndNormal;
        point       (0 0 0);
        normal      (0 0 1);
        interpolate yes;
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        planeType    | the method of specification of the plane | yes |
        interpolate  | interpolate values to the surface points | no | no
    \endtable

See also
    Foam::plane

SourceFiles
    sampledCutPlane.C

\*---------------------------------------------------------------------------*/

#ifndef sampledCutPlane_H
#define sampledCutPlane_H

#include "sampledIsoSurfaceSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSurfaces
{

/*---------------------------------------------------------------------------*\
                        Class cutPlane Declaration
\*---------------------------------------------------------------------------*/

class cutPlane
:
    public sampledIsoSurfaceSurface
{
    // Private Data

        //- Plane
        const plane plane_;


    // Private Member Functions

        //- Generate the iso surface
        virtual autoPtr<cutPolyIsoSurface> calcIsoSurf() const;


public:

    //- Runtime type information
    TypeName("cutPlane");


    // Constructors

        //- Construct from dictionary
        cutPlane
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~cutPlane();


    // Member Functions

        //- Does the surface need an update?
        virtual bool needsUpdate() const;

        //- Write
        virtual void print(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sampledSurfaces
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
