/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::coordinateSystems::coordinateSystems

Description
    Provides a centralised coordinateSystem collection.

    For example with the porous region specified in \c constant/fvModels as

    \verbatim
    porosity
    {
        type            explicitPorositySource;

        explicitPorositySourceCoeffs
        {
            select          cellZone;
            cellZone        porousBlockage;

            type            DarcyForchheimer;

            // D 100;  // Very little blockage
            // D 200;  // Some blockage but steady flow
            // D 500;  // Slight waviness in the far wake
            D 1000; // Fully shedding behavior

            d   ($D $D $D);
            f   (0 0 0);

            coordinateSystem porousBlockage;
        }
    }
    \endverbatim

    the corresponding coordinate system \c porousBlockage is looked-up
    automatically from the \c constant/coordinateSystems dictionary:

    \verbatim
    porousBlockage
    {
        type    cartesian;
        origin  (0 0 0);
        coordinateRotation
        {
            type    axesRotation;
            e1  (1 0 0);
            e2  (0 1 0);
        }
    }
    \endverbatim

    See \c tutorials/modules/incompressibleFluid/porousBlockage

SourceFiles
    coordinateSystems.C

\*---------------------------------------------------------------------------*/
#ifndef coordinateSystems_H
#define coordinateSystems_H

#include "coordinateSystem.H"
#include "PtrDictionary.H"
#include "DemandDrivenMeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coordinateSystems
{

/*---------------------------------------------------------------------------*\
                      Class coordinateSystems Declaration
\*---------------------------------------------------------------------------*/

class coordinateSystems
:
    public DemandDrivenMeshObject
    <
        objectRegistry,
        GeometricMeshObject,
        coordinateSystems
    >,
    public PtrDictionary<coordinateSystem>
{
protected:

    friend class DemandDrivenMeshObject
    <
        objectRegistry,
        GeometricMeshObject,
        coordinateSystems
    >;

    // Protected Constructors

        //- Read construct from objectRegistry
        explicit coordinateSystems(const objectRegistry& obr);


public:

    //- Runtime type information
    TypeName("coordinateSystems");


    // Constructors

        //- Disallow default bitwise copy construction
        coordinateSystems(const coordinateSystems&) = delete;

        using DemandDrivenMeshObject
        <
            objectRegistry,
            GeometricMeshObject,
            coordinateSystems
        >::New;


    //- Destructor
    virtual ~coordinateSystems()
    {}


    // Member Functions

        //- Is object global
        virtual bool global() const
        {
            return true;
        }

        //- ReadData function required for regIOobject read operation
        virtual bool readData(Istream&);

        //- WriteData function required for regIOobject write operation
        bool writeData(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const coordinateSystems&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coordinateSystems

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
struct typeGlobal
<
    DemandDrivenMeshObject
    <
        objectRegistry,
        GeometricMeshObject,
        coordinateSystems::coordinateSystems
    >
>
{
    static const bool global = true;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
