/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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
    Foam::epsilonLowReWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a turbulence dissipation wall function
    condition for low- and high-Reynolds number turbulent flow cases.

    The condition can be applied to wall boundaries, whereby it inserts near
    wall epsilon values directly into the epsilon equation to act as a
    constraint.

    The model operates in two modes, based on the computed laminar-to-turbulent
    switch-over y+ value derived from kappa and E.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        Cmu          | model coefficient       | no          | 0.09
        kappa        | Von Karman constant     | no          | 0.41
        E            | model coefficient       | no          | 9.8
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            epsilonLowReWallFunction;
    }
    \endverbatim

SeeAlso
    Foam::epsilonWallFunctionFvPatchScalarField

SourceFiles
    epsilonLowReWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef epsilonLowReWallFunctionFvPatchScalarField_H
#define epsilonLowReWallFunctionFvPatchScalarField_H

#include "epsilonWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class epsilonLowReWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class epsilonLowReWallFunctionFvPatchScalarField
:
    public epsilonWallFunctionFvPatchScalarField
{

protected:

    // Protected data

        //- Y+ at the edge of the laminar sublayer
        scalar yPlusLam_;


    // Protected Member Functions

        //- Calculate the Y+ at the edge of the laminar sublayer
        scalar yPlusLam(const scalar kappa, const scalar E);

        //- Calculate the epsilon and G
        virtual void calculate
        (
            const turbulenceModel& turbulence,
            const List<scalar>& cornerWeights,
            const fvPatch& patch,
            scalarField& G,
            scalarField& epsilon
        );


public:

    //- Runtime type information
    TypeName("epsilonLowReWallFunction");


    // Constructors

        //- Construct from patch and internal field
        epsilonLowReWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        epsilonLowReWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  epsilonLowReWallFunctionFvPatchScalarField
        //  onto a new patch
        epsilonLowReWallFunctionFvPatchScalarField
        (
            const epsilonLowReWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        epsilonLowReWallFunctionFvPatchScalarField
        (
            const epsilonLowReWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new epsilonLowReWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        epsilonLowReWallFunctionFvPatchScalarField
        (
            const epsilonLowReWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new epsilonLowReWallFunctionFvPatchScalarField(*this, iF)
            );
        }

    //- Destructor
    virtual ~epsilonLowReWallFunctionFvPatchScalarField()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
