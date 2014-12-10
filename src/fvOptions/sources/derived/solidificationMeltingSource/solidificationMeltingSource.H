/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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
    Foam::fv::solidificationMeltingSource

Description
    This source is designed to model the effect of solidification and melting
    processes, e.g. windhield defrosting.  The phase change occurs at the
    melting temperature, \c Tmelt.

    The presence of the solid phase in the flow field is incorporated into the
    model as a momentum porosity contribution; the energy associated with the
    phase change is added as an enthalpy contribution.

    Based on the references:

      1. V.R. Voller and C. Prakash, A fixed grid numerical modelling
         methodology for convection-diffusion mushy phase-change problems,
         Int. J. Heat Mass Transfer 30(8):17091719, 1987.
      2. C.R. Swaminathan. and V.R. Voller, A general enthalpy model for
         modeling solidification processes, Metallurgical Transactions
         23B:651664, 1992.

    The model generates a field \c \<name\>:alpha1 which can be visualised to
    to show the melt distribution as a fraction [0-1]

    \heading Source usage
    Example usage:
    \verbatim
    solidificationMeltingSource1
    {
        type            solidificationMeltingSource;
        active          on;
        selectionMode   cellZone;
        cellZone        iceZone;

        solidificationMeltingSourceCoeffs
        {
            Tmelt           273;
            L               334000;
            thermoMode      thermo;
            beta            50e-6;
            rhoRef          800;
        }
    }
    \endverbatim

    Where:
    \table
        Property     | Description             | Required    | Default value
        Tmelt        | Melting temperature [K] | yes         |
        L            | Latent heat of fusion [J/kg] | yes    |
        relax        | Relaxation coefficient [0-1] | no     | 0.9
        thermoMode   | Thermo mode [thermo|lookup] | yes     |
        rhoRef       | Reference (solid) density | yes       |
        rhoName      | Name of density field   | no          | rho
        TName        | Name of temperature field | no        | T
        CpName       | Name of specific heat capacity field | no | Cp
        UName        | Name of velocity field  | no          | U
        phiName      | Name of flux field      | no          | phi
        Cu           | Model coefficient       | no          | 100000
        q            | Model coefficient       | no          | 0.001
        beta         | Thermal expansion coefficient [1/K] | yes |
        g            | Accelerartion due to gravity | no     |
    \endtable

SourceFiles
    solidificationMeltingSource.C
    solidificationMeltingSourceIO.C

\*---------------------------------------------------------------------------*/

#ifndef solidificationMeltingSource_H
#define solidificationMeltingSource_H

#include "fvMesh.H"
#include "volFields.H"
#include "fvOption.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class solidificationMeltingSource Declaration
\*---------------------------------------------------------------------------*/

class solidificationMeltingSource
:
    public option
{
public:

    enum thermoMode
    {
        mdThermo,
        mdLookup
    };

    static const NamedEnum<thermoMode, 2> thermoModeTypeNames_;


private:

    // Private data

        //- Temperature at which melting occurs [K]
        scalar Tmelt_;

        //- Latent heat of fusion [J/kg]
        scalar L_;

        //- Phase fraction under-relaxation coefficient
        scalar relax_;

        //- Thermodynamics mode
        thermoMode mode_;

        //- Reference density - typically the solid density
        scalar rhoRef_;

        //- Name of temperature field - default = "T" (optional)
        word TName_;

        //- Name of specific heat capacity field - default = "Cp" (optional)
        word CpName_;

        //- Name of velocity field - default = "U" (optional)
        word UName_;

        //- Name of flux field - default = "phi" (optional)
        word phiName_;

        //- Mushy region momentum sink coefficient [1/s]; default = 10^5
        scalar Cu_;

        //- Coefficient used in porosity calc - default = 0.001
        scalar q_;

        //- Thermal expansion coefficient [1/K]
        scalar beta_;

        //- Phase fraction indicator field
        volScalarField alpha1_;

        //- Current time index (used for updating)
        label curTimeIndex_;

        //- Temperature change cached for source calculation when alpha1 updated
        scalarField deltaT_;


    // Private Member Functions

        //- Flag to indicate whether to solve for given field
        bool solveField(const word& fieldName) const;

        //- Return the specific heat capacity field
        tmp<volScalarField> Cp() const;

        //- Return the gravity vector
        vector g() const;

        //- Update the model
        void update(const volScalarField& Cp);

        //- Helper function to apply to the energy equation
        template<class RhoFieldType>
        void apply(const RhoFieldType& rho, fvMatrix<scalar>& eqn);

        //- Disallow default bitwise copy construct
        solidificationMeltingSource(const solidificationMeltingSource&);

        //- Disallow default bitwise assignment
        void operator=(const solidificationMeltingSource&);


public:

    //- Runtime type information
    TypeName("solidificationMeltingSource");


    // Constructors

        //- Construct from explicit source name and mesh
        solidificationMeltingSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        //- Flag to bypass the apply flag list checking
        virtual bool alwaysApply() const;


        // Add explicit and implicit contributions

            //- Add explicit contribution to enthalpy equation
            virtual void addSup(fvMatrix<scalar>& eqn, const label fieldI);

            //- Add implicit contribution to momentum equation
            virtual void addSup(fvMatrix<vector>& eqn, const label fieldI);


        // Add explicit and implicit contributions to compressible equation

            //- Add explicit contribution to compressible enthalpy equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldI
            );

            //- Add implicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );


        // I-O

            //- Write the source properties
            virtual void writeData(Ostream&) const;

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "solidificationMeltingSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
