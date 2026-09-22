/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::GurevichHSR

Description
    Heterogeneous Surface Reaction (HSR) model for aluminum particle combustion
    following Gurevich et al. (1970) and Zhang et al. (2022, Eqs. 10-15):

        Al + 0.75 O2 -> 0.5 Al2O3(condensed)

    - Effective reacting area: Ap,eff = pi * dp,Al^2
    - Enforces thermal ignition barrier: HSR activates once Tp >= Tmelt,Al (933 K)
    - Directly deposits newly formed condensed Al2O3 onto the particle
    - Depletes O2 from the surrounding Eulerian gas cell

\*---------------------------------------------------------------------------*/

#ifndef GurevichHSR_H
#define GurevichHSR_H

#include "SurfaceReactionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class GurevichHSR Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class GurevichHSR
:
    public SurfaceReactionModel<CloudType>
{
    // Private Data

        //- Pre-exponential rate coefficient [m/s]
        scalar Ar_;

        //- Activation energy [J/mol]
        scalar Ea_;

        //- Heat of surface reaction [J/kg Al]
        scalar h_HSR_;

        //- Aluminum melting temperature [K]
        scalar TmeltAl_;

        //- Density of solid/liquid aluminum [kg/m^3]
        scalar rhoAl_;


        // Addressing

        //- Al solid index in parcel
        mutable label AlSolidId_;

        //- Al2O3 solid index in parcel
        mutable label Al2O3SolidId_;

        //- Al liquid index in parcel
        mutable label AlLiquidId_;
        mutable label Al2O3LiquidId_;

        //- O2 carrier species index in Eulerian phase
        mutable label O2GlobalId_;


public:

    //- Runtime type information
    TypeName("GurevichHSR");


    // Constructors

        //- Construct from dictionary
        GurevichHSR
        (
            const dictionary& dict,
            CloudType& owner
        );

        //- Construct copy
        GurevichHSR
        (
            const GurevichHSR<CloudType>& srm
        );

        //- Construct and return a clone
        virtual autoPtr<SurfaceReactionModel<CloudType>> clone() const
        {
            return autoPtr<SurfaceReactionModel<CloudType>>
            (
                new GurevichHSR<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~GurevichHSR();


    // Member Functions

        //- Update surface reactions
        virtual scalar calculate
        (
            const scalar dt,
            const label celli,
            const scalar d,
            const scalar T,
            const scalar Tc,
            const scalar pc,
            const scalar rhoc,
            const scalar mass,
            const scalarField& YGas,
            const scalarField& YLiquid,
            const scalarField& YSolid,
            const scalarField& YMixture,
            const scalar N,
            scalarField& dMassGas,
            scalarField& dMassLiquid,
            scalarField& dMassSolid,
            scalarField& dMassSRCarrier
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "GurevichHSR.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
