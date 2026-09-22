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

\*---------------------------------------------------------------------------*/

#include "OxideCapEvaporation.H"
#include "specie.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::tmp<Foam::scalarField> Foam::OxideCapEvaporation<CloudType>::calcXc
(
    const label celli
) const
{
    scalarField Xc(this->owner().composition().carrier().Y().size());

    forAll(Xc, i)
    {
        Xc[i] =
            this->owner().composition().carrier().Y()[i][celli]
           /this->owner().composition().carrier().WiValue(i);
    }

    return Xc/max(sum(Xc), scalar(1e-12));
}


template<class CloudType>
Foam::scalar Foam::OxideCapEvaporation<CloudType>::Sh
(
    const scalar Re,
    const scalar Sc
) const
{
    return 2.0 + 0.6*Foam::sqrt(max(Re, scalar(0.0)))*cbrt(max(Sc, scalar(1e-4)));
}


template<class CloudType>
Foam::scalar Foam::OxideCapEvaporation<CloudType>::calcBeta
(
    const scalar Tp
) const
{
    // If dictionary explicitly specifies a constant beta, honor it
    if (this->coeffDict().found("beta"))
    {
        return this->coeffDict().template lookup<scalar>("beta");
    }

    // Before alumina melts (Tp < 2327 K), oxide is an intact shell hindering bulk evaporation
    if (Tp < TmeltOxide_)
    {
        return 0.95; // Heavily restricted until oxide shell melts into a cap
    }

    // Above 2327 K: Oxide melts into a cap on top of liquid aluminum core
    // beta = 0.5 * (1 - cos(alpha)) (Eq. 32)
    return 0.20;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::OxideCapEvaporation<CloudType>::OxideCapEvaporation
(
    const dictionary& dict,
    CloudType& owner
)
:
    PhaseChangeModel<CloudType>(dict, owner, typeName),
    liquids_(owner.thermo().liquids()),
    condensation_
    (
        this->coeffDict().template lookupOrDefault<Switch>
        (
            "condensation",
            false
        )
    ),
    activeLiquids_(this->coeffDict().lookup("activeLiquids")),
    liqToCarrierMap_(activeLiquids_.size(), -1),
    liqToLiqMap_(activeLiquids_.size(), -1),
    A_clausius_(this->coeffDict().template lookupOrDefault<scalar>("A_clausius", 12.19)),
    B_clausius_(this->coeffDict().template lookupOrDefault<scalar>("B_clausius", 34037.0)),
    Pref_(this->coeffDict().template lookupOrDefault<scalar>("Pref", 101325.0)),
    TmeltOxide_(this->coeffDict().template lookupOrDefault<scalar>("TmeltOxide", 2327.0)),
    L_evap_(this->coeffDict().template lookupOrDefault<scalar>("L_evap", 1.05e7))
{
    if (activeLiquids_.size() == 0)
    {
        WarningInFunction
            << "Evaporation model selected, but no active liquids defined"
            << nl << endl;
    }
    else
    {
        Info<< "Participating liquid species:" << endl;

        forAll(activeLiquids_, i)
        {
            Info<< "    " << activeLiquids_[i] << endl;
            if (owner.composition().carrier().species().found(activeLiquids_[i]))
            {
                liqToCarrierMap_[i] =
                    owner.composition().carrierId(activeLiquids_[i]);
            }

            const label idLiquid = owner.composition().idLiquid();
            if (idLiquid >= 0)
            {
                liqToLiqMap_[i] =
                    owner.composition().localId(idLiquid, activeLiquids_[i]);
            }
        }
    }
}


template<class CloudType>
Foam::OxideCapEvaporation<CloudType>::OxideCapEvaporation
(
    const OxideCapEvaporation<CloudType>& pcm
)
:
    PhaseChangeModel<CloudType>(pcm),
    liquids_(pcm.owner().thermo().liquids()),
    condensation_(pcm.condensation_),
    activeLiquids_(pcm.activeLiquids_),
    liqToCarrierMap_(pcm.liqToCarrierMap_),
    liqToLiqMap_(pcm.liqToLiqMap_),
    A_clausius_(pcm.A_clausius_),
    B_clausius_(pcm.B_clausius_),
    Pref_(pcm.Pref_),
    TmeltOxide_(pcm.TmeltOxide_),
    L_evap_(pcm.L_evap_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::OxideCapEvaporation<CloudType>::~OxideCapEvaporation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::OxideCapEvaporation<CloudType>::calculate
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar Re,
    const scalar Pr,
    const scalar d,
    const scalar nu,
    const scalar T,
    const scalar Ts,
    const scalar pc,
    const scalar Tc,
    const scalarField& X,
    scalarField& dMassPC
) const
{
    // Construct carrier phase species mole fractions
    const scalarField Xc(calcXc(p.cell()));

    forAll(activeLiquids_, i)
    {
        const label gid = liqToCarrierMap_[i];
        const label lid = liqToLiqMap_[i];

        if (lid < 0 || gid < 0) continue;

        // 1. Vapor Pressure from Clausius-Clapeyron / Alcock formulation (Eq. 19)
        const scalar exponent = A_clausius_ - B_clausius_ / max(T, scalar(300.0));
        scalar pSat = Pref_ * Foam::exp(min(exponent, scalar(15.0)));

        // Vapor mole fraction at droplet surface from Raoult's law (Eq. 19)
        const scalar XFs = min(pSat / max(pc, scalar(1e4)), scalar(0.99));

        // Molecular weight of fuel vapor and non-vapor carrier gas
        const scalar WF = liquids_.properties()[lid].W();
        const scalar WnonF = 28.84; // Average molecular weight of carrier (Air)

        // Vapor mass fraction at droplet surface (Eq. 18)
        const scalar YFs = (XFs * WF) / max(XFs * WF + (1.0 - XFs) * WnonF, scalar(1e-6));

        // Far-field vapor mass fraction in Eulerian cell
        const scalar YFinf = max(this->owner().composition().carrier().Y()[gid][p.cell()], scalar(0.0));

        // Spalding mass transfer number (Eq. 17)
        const scalar BM = max((YFs - YFinf) / max(1.0 - YFs, scalar(1e-4)), scalar(0.0));

        if (BM <= 1e-12) continue;

        // Vapor diffusivity [m^2/s]
        // If liquid properties specify unphysical liquid-state D (< 1e-6),
        // fallback to gas-phase binary diffusion with Sc ~ 0.7
        scalar Dab = liquids_.properties()[lid].D(pc, Ts);
        if (Dab < 1e-6)
        {
            Dab = (nu + rootVSmall) / 0.7;
        }
        const scalar Sc = nu / (Dab + rootVSmall);
        const scalar Sh = this->Sh(Re, Sc);

        // Core effective diameter
        const scalar dpAl = d;

        // 2. Oxide cap blocking fraction beta (Eq. 32)
        const scalar beta = calcBeta(T);

        // 3. Evaporation rate: mDot = pi * (1 - beta) * dpAl * Sh * rho * DF * ln(1 + BM) (Eq. 16)
        const scalar rhoc = td.rhoc();
        const scalar mDotEvap = pi * (1.0 - beta) * dpAl * Sh * rhoc * Dab * Foam::log(1.0 + BM);

        scalar dMass = mDotEvap * dt;
        if (!condensation_)
        {
            dMass = max(dMass, scalar(0.0));
        }

        // Add to parcel phase change transfer
        dMassPC[lid] += dMass;
    }
}


template<class CloudType>
Foam::scalar Foam::OxideCapEvaporation<CloudType>::dh
(
    const label idc,
    const label idl,
    const scalar p,
    const scalar T
) const
{
    return L_evap_;
}


template<class CloudType>
Foam::scalar Foam::OxideCapEvaporation<CloudType>::Tvap
(
    const scalarField& X
) const
{
    return liquids_.Tpt(X);
}


template<class CloudType>
Foam::scalar Foam::OxideCapEvaporation<CloudType>::TMax
(
    const scalar p,
    const scalarField& X
) const
{
    return 2792.0; // Aluminum normal boiling point
}

// ************************************************************************* //
