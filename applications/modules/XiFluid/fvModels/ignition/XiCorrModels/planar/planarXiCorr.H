/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::XiCorrModels::planar

Description
    Planar ignition kernel flame-wrinkling correction model

Usage
    Example usage:
    \verbatim
    XiCorr
    {
        type            planar;
        select          all;
        area            1e-4;
    }
    \endverbatim

    Where:
    \table
        Property   | Description                   | Required   | Default value
        select     | Correction cell set selection | yes        |
        area       | Cross-sectional area of the 1D mesh | yes  |
        bMin       | Min b below which no correction | no | 0.01
    \endtable

SourceFiles
    planar.C

\*---------------------------------------------------------------------------*/

#ifndef planarXiCorr_H
#define planarXiCorr_H

#include "XiCorrModel.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiCorrModels
{

/*---------------------------------------------------------------------------*\
                           Class planar Declaration
\*---------------------------------------------------------------------------*/

class planar
:
    public XiCorrModel
{
    // Private Data

        //- Cross-sectional area of the 1D mesh
        dimensionedScalar area_;

        //- Update coefficients from given dictionary
        virtual bool readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("planar");


    // Constructors

        //- Construct from mesh and dictionary
        planar(const fvMesh& mesh, const dictionary& dict);

        //- Disallow default bitwise copy construction
        planar(const planar&) = delete;


    //- Destructor
    virtual ~planar();


    // Member Functions

        //- Return the area of the ignition kernel
        //  calculated from the given kernel volume
        virtual dimensionedScalar Ak(const dimensionedScalar& Vk) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const planar&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiCorrModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
