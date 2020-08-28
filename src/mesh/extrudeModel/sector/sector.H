/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::extrudeModels::sector

Description
    Extrudes by rotating a surface around an axis
    - extrusion is opposite the surface/patch normal so inwards the source
      mesh
    - axis direction has to be consistent with this.
    - use -mergeFaces option if doing full 360 and want to merge front and back
    - note direction of axis. This should be consistent with rotating against
      the patch normal direction. If you get it wrong you'll see all cells
      with extreme aspect ratio and internal faces wrong way around in
      checkMesh

\*---------------------------------------------------------------------------*/

#ifndef sector_H
#define sector_H

#include "extrudeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace extrudeModels
{

/*---------------------------------------------------------------------------*\
                           Class sector Declaration
\*---------------------------------------------------------------------------*/

class sector
:
    public extrudeModel
{
    // Private Data

        //- Point on axis
        const point axisPt_;

        //- Normalized direction of axis
        const vector axis_;

        //- Overall angle (radians)
        const scalar angle_;


public:

    //- Runtime type information
    TypeName("sector");

    // Constructors

        //- Construct from dictionary
        sector(const dictionary& dict);


    //- Destructor
    virtual ~sector();


    // Member Operators

        point operator()
        (
            const point& surfacePoint,
            const vector& surfaceNormal,
            const label layer
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace extrudeModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
