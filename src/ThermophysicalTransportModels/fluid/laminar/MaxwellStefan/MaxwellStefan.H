/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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
    Foam::turbulenceThermophysicalTransportModels::MaxwellStefan

Description
    Base class for multi-component Maxwell Stefan generalised Fick's law
    diffusion coefficients based temperature gradient heat flux model with
    optional Soret thermal diffusion of species.

    The binary mass diffusion coefficients are specified as Function2<scalar>s
    of pressure and temperature but independent of composition.

    The heat flux source is implemented as an implicit energy correction to the
    temperature gradient based flux source.  At convergence the energy
    correction is 0.

    References:
    \verbatim
        Taylor, R., & Krishna, R. (1993).
        Multicomponent mass transfer (Vol. 2).
        John Wiley & Sons.

        Merk, H. J. (1959).
        The macroscopic equations for simultaneous heat and mass transfer
        in isotropic, continuous and closed systems.
        Applied Scientific Research,
        Section A, 8(1), 73-99.
    \endverbatim

SourceFiles
    MaxwellStefan.C

\*---------------------------------------------------------------------------*/

#include "Function2.H"
#include "LUscalarMatrix.H"

#ifndef MaxwellStefan_H
#define MaxwellStefan_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class MaxwellStefan Declaration
\*---------------------------------------------------------------------------*/

template<class BasicThermophysicalTransportModel>
class MaxwellStefan
:
    public BasicThermophysicalTransportModel,
    public UpdateableMeshObject<fvMesh>
{
    // Private data

        // Model coefficients

            //- Array of specie binary mass diffusion coefficient functions
            //  [m^2/s]
            List<PtrList<Function2<scalar>>> DFuncs_;

            //- List of specie Soret thermal diffusion coefficient
            //  functions [kg/m/s]
            PtrList<Function2<scalar>> DTFuncs_;

            //- Generalised Fick's law diffusion coefficients field list.
            //  This is the diagonal of the mass diffusion coefficient matrix
            mutable PtrList<volScalarField> Dii_;

            //- List of fields of the explicit part of the mass diffusion flux
            // of the species
            mutable PtrList<surfaceScalarField> jexp_;


        // Workspace for diffusion coefficient transformation

            //- Molecular weights of the species
            mutable scalarField W;

            //- List of mass-fraction field pointers
            //  for the internal mesh or a patch
            mutable List<const scalarField*> YPtrs;

            //- Matrix of mass diffusion coefficient field pointers
            //  for the internal mesh or a patch
            mutable SquareMatrix<scalarField*> DijPtrs;

            //- Mass-fractions at a cell or face
            mutable scalarField Y;

            //- Mole-fractions at a cell or face
            mutable scalarField X;

            //- Binary mass diffusion coefficient matrix at a cell or face
            mutable scalarSquareMatrix DD;

            //- Matrix form of the coefficients in the Maxwell-Stefan equation
            //  at a cell or face
            mutable LUscalarMatrix A;

            //- Matrix of composition coefficients at a cell or face
            //  used to transform the binary mass diffusion coefficients into
            //  the generalised Fick's law diffusion coefficients
            mutable scalarSquareMatrix B;

            //- Inverse of A
            mutable scalarSquareMatrix invA;

            //- Matrix of the generalised Fick's law diffusion coefficients
            //  at a cell or face
            mutable scalarSquareMatrix D;


    // Private member functions

        //- Transform the Binary mass diffusion coefficient matrix DD into the
        //  matrix of the generalised Fick's law diffusion coefficients D
        void transformDiffusionCoefficient() const;

        //- Calls transformDiffusionCoefficient() for each cell and patch face
        //  and maps the resulting D into DijPtrs
        void transformDiffusionCoefficientFields() const;

        //- Calls transformDiffusionCoefficientFields() for the internal mesh
        //  and each patch and maps the resulting DijPtrs into Dii_ and Dij
        void transform(List<PtrList<volScalarField>>& Dij) const;

        //- Update Dii_
        void updateDii() const;

        //- Update Dii_ if not yet constructed and return
        const PtrList<volScalarField>& Dii() const;

        //- Update jexp_ if not yet constructed and return
        const PtrList<surfaceScalarField>& jexp() const;


public:

    typedef typename BasicThermophysicalTransportModel::alphaField
        alphaField;

    typedef typename
        BasicThermophysicalTransportModel::momentumTransportModel
        momentumTransportModel;

    typedef typename BasicThermophysicalTransportModel::thermoModel
        thermoModel;


    // Constructors

        //- Construct from a momentum transport model and a thermo model
        MaxwellStefan
        (
            const word& type,
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );


    //- Destructor
    virtual ~MaxwellStefan()
    {}


    // Member Functions

        //- Read thermophysicalTransport dictionary
        virtual bool read();

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction [kg/m/s]
        //  This is the self-diffusion coefficient
        virtual tmp<volScalarField> DEff(const volScalarField& Yi) const;

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction for patch [kg/m/s]
        virtual tmp<scalarField> DEff
        (
            const volScalarField& Yi,
            const label patchi
        ) const;

        //- Return the heat flux [W/m^2]
        virtual tmp<surfaceScalarField> q() const;

        //- Return the source term for the energy equation
        virtual tmp<fvScalarMatrix> divq(volScalarField& he) const;

        //- Return the specie flux for the given specie mass-fraction [kg/m^2/s]
        virtual tmp<surfaceScalarField> j(const volScalarField& Yi) const;

        //- Return the source term for the given specie mass-fraction equation
        virtual tmp<fvScalarMatrix> divj(volScalarField& Yi) const;

        //- Update the diffusion coefficients and flux corrections
        virtual void predict();


        // Mesh changes

            //- Update for mesh motion
            virtual bool movePoints();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap& map);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap& map);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap& map);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MaxwellStefan.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
