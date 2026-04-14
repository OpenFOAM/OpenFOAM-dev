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
    Foam::zonalThermo

Description
    Base-class for zonal thermodynamic properties.

See also
    Foam::basicThermo

SourceFiles
    zonalThermo.C

\*---------------------------------------------------------------------------*/

#ifndef zonalThermo_H
#define zonalThermo_H

#include "basicThermo.H"
#include "zonalThermoZones.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class zonalThermo Declaration
\*---------------------------------------------------------------------------*/

class zonalThermo
:
    virtual public basicThermo
{
public:

    //- Runtime type information
    TypeName("zonalThermo");


    // Public Classes

        //- Forward declare the implementation class
        class implementation;


    //- Destructor
    virtual ~zonalThermo();
};


/*---------------------------------------------------------------------------*\
                 Class zonalThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class zonalThermo::implementation
:
    virtual public zonalThermo
{
protected:

    // Protected data

        //- Addressing
        const zonalThermoZones& ztzs_;


public:

    // Constructors

        //- Construct from dictionary, zone names, mesh and phase name
        implementation
        (
            const dictionary&,
            const hashedWordList&,
            const fvMesh&,
            const word&
        );


    //- Destructor
    virtual ~implementation();


    // Member Functions

        //- Get the slicer
        inline nil Yslicer() const;

        //- Get the composition of an internal cell
        inline label cellComposition
        (
            const nil,
            const label celli
        ) const;

        //- Get the composition of a boundary face
        inline label patchFaceComposition
        (
            const nil,
            const label patchi,
            const label patchFacei
        ) const;

        //- Get the slicer for the given source
        inline nil Yslicer
        (
            const fvSource& model,
            const volScalarField::Internal& source
        ) const;

        //- Get the composition of a source cell
        inline label sourceCellComposition
        (
            const nil,
            const label i
        ) const;

        //- Get the slicer for the given source
        inline const labelUList& Yslicer
        (
            const fvSource& model,
            const scalarField& source,
            const labelUList& cells
        ) const;

        //- Get the composition of a source cell
        inline label sourceCellComposition
        (
            const labelUList& cells,
            const label i
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "zonalThermoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
