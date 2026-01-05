/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::Lagrangian::RanzMarshallHeatTransfer

Description
    Heat transfer model of Ranz and Marshall for spherical particles in a
    turbulent carrier flow

Usage
    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type        RanzMarshallHeatTransfer;
    }
    \endverbatim

SourceFiles
    RanzMarshallHeatTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef RanzMarshallHeatTransfer_H
#define RanzMarshallHeatTransfer_H

#include "heatTransfer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{

/*---------------------------------------------------------------------------*\
                          Class RanzMarshallHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class RanzMarshallHeatTransfer
:
    public heatTransfer
{
private:

    // Private Member Functions

        //- Calculate the heat transfer coefficient
        virtual tmp<LagrangianSubScalarField> calcH
        (
            const LagrangianModelRef& model,
            const LagrangianSubMesh& subMesh
        ) const;


public:

    //- Runtime type information
    TypeName("RanzMarshallHeatTransfer");


    // Constructors

        //- Inherit constructors
        using heatTransfer::heatTransfer;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Lagrangian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
