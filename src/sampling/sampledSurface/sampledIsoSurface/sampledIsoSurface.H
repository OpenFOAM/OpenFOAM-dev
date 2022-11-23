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
    Foam::sampledSurfaces::isoSurface

Description
    A sampledSurface defined by an iso-value of a field.

    Example:
    \verbatim
    {
        type        isoSurface;
        isoField    p;
        isoValue    1e5;
        interpolate yes;
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        isoField     | the field to get an iso-surface of | yes |
        isoValue     | the iso-surface value   | yes         |
        interpolate  | interpolate values to the surface points | no | no
    \endtable

SourceFiles
    sampledIsoSurface.C

\*---------------------------------------------------------------------------*/

#ifndef sampledIsoSurface_H
#define sampledIsoSurface_H

#include "sampledIsoSurfaceSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSurfaces
{

/*---------------------------------------------------------------------------*\
                         Class isoSurface Declaration
\*---------------------------------------------------------------------------*/

class isoSurface
:
    public sampledIsoSurfaceSurface
{
    // Private Data

        //- Field to get iso surface of
        const word isoFieldName_;

        //- Iso values
        const scalarField isoValues_;


    // Private Member Functions

        //- Generate the iso surface
        virtual autoPtr<cutPolyIsoSurface> calcIsoSurf() const;


public:

    //- Runtime type information
    TypeName("isoSurface");


    // Constructors

        //- Construct from dictionary
        isoSurface
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~isoSurface();


    // Member Functions

        //- Return the list of fields required
        virtual wordList fields() const;

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
