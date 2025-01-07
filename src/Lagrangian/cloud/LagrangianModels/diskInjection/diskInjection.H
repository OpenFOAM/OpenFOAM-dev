/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::diskInjection

Description
    Disk injection model. This injects particles continuously over a disk with
    a given number rate. The disk is characterised by a centre, and axis and an
    inner and an outer diameter. The centre and axis and number rate are all
    Function1-s and can vary with time.

    Note that this model only controls the number and position of injected
    Lagrangian particles. All physical properties are specified by
    corresponding source conditions. So the velocity/direction/angle/etc..., is
    controlled by the velocity source condition, the size distribution by the
    diameter source condition, and the flow rate by the number source
    condition.

Usage
    \table
        Property      | Description                       | Required? | Default
        centre        | The centre of the disk            | yes       |
        axis          | The axis normal to the disk       | yes       |
        diameter      | The diameter of the circle        | if innerDiameter \
                                                            and outerDiameter \
                                                            are not specified |
        innerDiameter | The inner diameter of the annulus | if diameter is not \
                                                            specified |
        outerDiameter | The outer diameter of the annulus | if diameter is not \
                                                            specified |
        numberRate    | The number of particles to \
                        inject per unit time              | yes       |
    \endtable

    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type        pointInjection;
        point       (1 2 3) [cm];
        numberRate  10000;
    }
    \endverbatim

See also
    Foam::Function1s

SourceFiles
    diskInjection.C

\*---------------------------------------------------------------------------*/

#ifndef diskInjection_H
#define diskInjection_H

#include "LagrangianInjection.H"
#include "Function1.H"
#include "restartableRandomGenerator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{

/*---------------------------------------------------------------------------*\
                          Class diskInjection Declaration
\*---------------------------------------------------------------------------*/

class diskInjection
:
    public LagrangianInjection
{
private:

    // Private Data

        //- The centre of the disk
        autoPtr<Function1<point>> centre_;

        //- The axis of the disk
        autoPtr<Function1<vector>> axis_;

        //- The inner diameter
        scalar innerDiameter_;

        //- The outer diameter
        scalar outerDiameter_;

        //- The rate at which to inject
        autoPtr<Function1<scalar>> numberRate_;

        //- The number deferred to the next injection step
        mutable scalar numberDeferred_;

        //- A random number generator
        mutable restartableRandomGenerator rndGen_;

        //- The time index
        mutable label timeIndex_;

        //- Cached radius fractions
        mutable autoPtr<scalarField> rFracPtr_;

        //- Cached axes
        mutable autoPtr<vectorField> axisPtr_;

        //- Cached radial vectors
        mutable autoPtr<vectorField> radialPtr_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& modelDict);


public:

    //- Runtime type information
    TypeName("diskInjection");


    // Constructors

        //- Construct from components
        diskInjection
        (
            const word& name,
            const LagrangianMesh& mesh,
            const dictionary& modelDict,
            const dictionary& stateDict
        );


    // Member Functions

        // Access

            //- Access the cached radius fractions. Only valid during injection.
            const scalarField& rFrac() const;

            //- Access the cached axes. Only valid during injection.
            const vectorField& axis() const;

            //- Access the cached radial vectors. Only valid during injection.
            const vectorField& radial() const;


        //- Create new elements in the Lagrangian mesh
        LagrangianSubMesh modify
        (
            LagrangianMesh& mesh,
            const LagrangianSubMesh&
        ) const;


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& modelDict);

            //- Write state
            virtual void writeState(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Lagrangian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
