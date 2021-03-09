/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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
    Foam::fv::VoFSolidificationMeltingSource

Description
    Solidification and melting model for VoF simulations.

    The presence of the solid phase in the flow field is incorporated into the
    model as a momentum porosity contribution; the energy associated with the
    phase change is added as an enthalpy contribution.  The solid fraction as a
    function of temperature \c alphaSolidT is specified as a Foam::Function1.

    The model writes the field \c alpha[01].solid which can be visualised to to
    show the solid distribution.

Usage
    Example usage:
    \verbatim
    VoFSolidificationMeltingSource1
    {
        type            VoFSolidificationMeltingSource;

        selectionMode   cellZone;
        cellZone        solidZone;

        alphaSolidT     table
        (
            (330 1)
            (335 0)
        );

        L               334000;
    }
    \endverbatim

    Where:
    \table
        Property     | Description             | Required    | Default value
        alphaSolidT  | Solid fraction as function of temperature | yes |
        L            | Latent heat of fusion [J/kg] | yes    |
        relax        | Relaxation coefficient [0-1] | no     | 0.9
        Cu           | Model coefficient       | no          | 100000
        q            | Model coefficient       | no          | 0.001
    \endtable

See also
    Foam::fv::solidificationMeltingSource
    Foam::Function1

SourceFiles
    VoFSolidificationMeltingSource.C
    VoFSolidificationMeltingSourceIO.C

\*---------------------------------------------------------------------------*/

#ifndef VoFSolidificationMeltingSource_H
#define VoFSolidificationMeltingSource_H

#include "fvModel.H"
#include "fvCellSet.H"
#include "fvMesh.H"
#include "volFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class VoFSolidificationMeltingSource Declaration
\*---------------------------------------------------------------------------*/

class VoFSolidificationMeltingSource
:
    public fvModel
{
    // Private Data

        //- The set of cells the fvConstraint applies to
        fvCellSet set_;

        //- Solid fraction as a function of temperature
        autoPtr<Function1<scalar>> alphaSolidT_;

        //- Latent heat of fusion [J/kg]
        dimensionedScalar L_;

        //- Phase fraction under-relaxation coefficient
        scalar relax_;

        //- Mushy region momentum sink coefficient [1/s]; default = 10^5
        scalar Cu_;

        //- Coefficient used in porosity calc - default = 0.001
        scalar q_;

        //- Solid phase fraction
        mutable volScalarField alphaSolid_;

        //- Current time index (used for updating)
        mutable label curTimeIndex_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();

        //- Return the name of the solid phase fraction
        word alphaSolidName() const;

        //- Update the model
        void update() const;

        //- Helper function to apply to the energy equation
        template<class RhoFieldType>
        void apply(const RhoFieldType& rho, fvMatrix<scalar>& eqn) const;


public:

    //- Runtime type information
    TypeName("VoFSolidificationMeltingSource");


    // Constructors

        //- Construct from explicit source name and mesh
        VoFSolidificationMeltingSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        VoFSolidificationMeltingSource
        (
            const VoFSolidificationMeltingSource&
        ) = delete;


    // Member Functions

        // Checks

            //- Return the list of fields for which the option adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Add explicit and implicit contributions to compressible equation

            //- Add explicit contribution to compressible enthalpy equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const word& fieldName
            ) const;

            //- Add implicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const word& fieldName
            ) const;


        // Mesh motion

            //- Update for mesh changes
            virtual void updateMesh(const mapPolyMesh&);


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const VoFSolidificationMeltingSource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
