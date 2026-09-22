/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "ReactingMultiphaseParcel.H"
#include "CompositionModel.H"
#include "NoDevolatilisation.H"
#include "NoSurfaceReaction.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;


// * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::CpEff
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[idG]*cloud.composition().Cp(idG, YGas_, p, T)
      + this->Y_[idL]*cloud.composition().Cp(idL, YLiquid_, p, T)
      + this->Y_[idS]*cloud.composition().Cp(idS, YSolid_, p, T);
}


template<class ParcelType>
template<class TrackCloudType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::hsEff
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[idG]*cloud.composition().hs(idG, YGas_, p, T)
      + this->Y_[idL]*cloud.composition().hs(idL, YLiquid_, p, T)
      + this->Y_[idS]*cloud.composition().hs(idS, YSolid_, p, T);
}


template<class ParcelType>
template<class TrackCloudType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::LEff
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[idG]*cloud.composition().L(idG, YGas_, p, T)
      + this->Y_[idL]*cloud.composition().L(idL, YLiquid_, p, T)
      + this->Y_[idS]*cloud.composition().L(idS, YSolid_, p, T);
}


template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::updateMassFractions
(
    const scalar mass0,
    const scalarField& dMassGas,
    const scalarField& dMassLiquid,
    const scalarField& dMassSolid,
    const label idG,
    const label idL,
    const label idS
)
{
    scalarField& YMix = this->Y_;

    scalar massGas =
        this->updateMassFraction(mass0*YMix[idG], dMassGas, YGas_);
    scalar massLiquid =
        this->updateMassFraction(mass0*YMix[idL], dMassLiquid, YLiquid_);
    scalar massSolid =
        this->updateMassFraction(mass0*YMix[idS], dMassSolid, YSolid_);

    scalar massNew = max(massGas + massLiquid + massSolid, rootVSmall);

    YMix[idG] = massGas/massNew;
    YMix[idL] = massLiquid/massNew;
    YMix[idS] = 1.0 - YMix[idG] - YMix[idL];

    return massNew;
}


// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::setCellValues(cloud, td);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::cellValueSourceCorrection
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    // Reuse correction from reacting parcel
    ParcelType::cellValueSourceCorrection(cloud, td, dt);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    typedef typename TrackCloudType::thermoCloudType thermoCloudType;
    const CompositionModel<thermoCloudType>& composition =
        cloud.composition();


    // Define local properties at beginning of timestep
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const scalar np0 = this->nParticle_;
    const scalar d0 = this->d_;
    const vector& U0 = this->U_;
    const scalar T0 = this->T_;
    const scalar mass0 = this->mass();

    const scalar pc = td.pc();

    const scalarField& YMix = this->Y_;
    const label idG = composition.idGas();
    const label idL = composition.idLiquid();
    const label idS = composition.idSolid();


    // Calc surface values
    scalar Ts, rhos, mus, Prs, kappas;
    this->calcSurfaceValues(cloud, td, T0, Ts, rhos, mus, Prs, kappas);
    scalar Res = this->Re(rhos, U0, td.Uc(), d0, mus);


    // Sources
    //~~~~~~~~

    // Explicit momentum source for particle
    vector Su = Zero;

    // Linearised momentum source coefficient
    scalar Spu = 0.0;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = Zero;

    // Explicit enthalpy source for particle
    scalar Sh = 0.0;

    // Linearised enthalpy source coefficient
    scalar Sph = 0.0;

    // Sensible enthalpy transfer from the particle to the carrier phase
    scalar dhsTrans = 0.0;


    // 1. Compute models that contribute to mass transfer - U, T held constant
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Phase change in liquid phase
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Mass transfer due to phase change
    scalarField dMassPC(YLiquid_.size(), 0.0);

    // Molar flux of species emitted from the particle (kmol/m^2/s)
    scalar Ne = 0.0;

    // Sum Ni*Cpi*Wi of emission species
    scalar NCpW = 0.0;

    // Surface concentrations of emitted species
    scalarField Cs(composition.carrier().species().size(), 0.0);

    // Calc mass and enthalpy transfer due to phase change
    this->calcPhaseChange
    (
        cloud,
        td,
        dt,
        Res,
        Prs,
        Ts,
        mus/rhos,
        d0,
        T0,
        mass0,
        idL,
        YMix[idL],
        YLiquid_,
        dMassPC,
        Sh,
        Ne,
        NCpW,
        Cs
    );


    // Devolatilisation
    // ~~~~~~~~~~~~~~~~

    // Mass transfer due to devolatilisation
    scalarField dMassDV(YGas_.size(), 0.0);

    // Calc mass and enthalpy transfer due to devolatilisation
    calcDevolatilisation
    (
        cloud,
        td,
        dt,
        Ts,
        d0,
        T0,
        mass0,
        mass0_,
        YMix[idG]*YGas_,
        YMix[idL]*YLiquid_,
        YMix[idS]*YSolid_,
        canCombust_,
        dMassDV,
        Sh,
        Ne,
        NCpW,
        Cs
    );


    // Surface reactions
    // ~~~~~~~~~~~~~~~~~

    // Change in carrier phase composition due to surface reactions
    scalarField dMassSRGas(YGas_.size(), 0.0);
    scalarField dMassSRLiquid(YLiquid_.size(), 0.0);
    scalarField dMassSRSolid(YSolid_.size(), 0.0);
    scalarField dMassSRCarrier(composition.carrier().species().size(), 0.0);

    // Calc mass and enthalpy transfer due to surface reactions
    calcSurfaceReactions
    (
        cloud,
        td,
        dt,
        d0,
        T0,
        mass0,
        canCombust_,
        Ne,
        YMix,
        YGas_,
        YLiquid_,
        YSolid_,
        dMassSRGas,
        dMassSRLiquid,
        dMassSRSolid,
        dMassSRCarrier,
        Sh,
        dhsTrans
    );


    // 2. Update the parcel properties due to change in mass
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Physical constraint: Total liquid consumption in dt cannot exceed available liquid mass
    if (idL >= 0)
    {
        const scalar mLiqAvail = mass0 * YMix[idL];
        const scalar totLiqLoss = sum(dMassPC) + sum(dMassSRLiquid);
        if (totLiqLoss > mLiqAvail && mLiqAvail > 0.0)
        {
            const scalar scale = mLiqAvail / max(totLiqLoss, rootVSmall);
            dMassPC *= scale;
            dMassSRLiquid *= scale;
        }
        else if (mLiqAvail <= 0.0)
        {
            dMassPC = 0.0;
            dMassSRLiquid = 0.0;
        }
    }

    scalarField dMassGas(dMassDV + dMassSRGas);
    scalarField dMassLiquid(dMassPC + dMassSRLiquid);
    scalarField dMassSolid(dMassSRSolid);
    scalar mass1 =
        updateMassFractions
        (
            mass0,
            dMassGas,
            dMassLiquid,
            dMassSolid,
            idG,
            idL,
            idS
        );

    this->Cp_ = CpEff(cloud, td, pc, T0, idG, idL, idS);

    // Update particle density or diameter
    if (cloud.constProps().constantVolume())
    {
        this->rho_ = mass1/this->volume();
    }
    else
    {
        this->d_ = cbrt(mass1/this->rho_*6.0/pi);
    }

    // Remove the particle when mass falls below minimum threshold
    if (np0*mass1 < cloud.constProps().minParcelMass())
    {
        td.keepParticle = false;

        if (cloud.solution().coupled())
        {
            scalar dm = np0*mass0;

            // Absorb parcel into carrier phase
            forAll(YGas_, i)
            {
                label gid = composition.localToCarrierId(idG, i);
                cloud.rhoTrans(gid)[this->cell()] += dm*YMix[idG]*YGas_[i];
            }
            forAll(YLiquid_, i)
            {
                label gid = composition.localToCarrierId(idL, i);
                cloud.rhoTrans(gid)[this->cell()] += dm*YMix[idL]*YLiquid_[i];
            }

            // No mapping between solid components and carrier phase
            /*
            forAll(YSolid_, i)
            {
                label gid = composition.localToCarrierId(idS, i);
                cloud.rhoTrans(gid)[this->cell()] += dm*YMix[idS]*YSolid_[i];
            }
            */

            cloud.UTransRef()[this->cell()] += dm*U0;

            cloud.hsTransRef()[this->cell()] +=
                dm*hsEff(cloud, td, pc, T0, idG, idL, idS);

            cloud.phaseChange().addToPhaseChangeMass(np0*mass1);
        }

        return;
    }

    // Correct surface values due to emitted species
    this->correctSurfaceValues(cloud, td, Ts, Cs, rhos, mus, Prs, kappas);
    Res = this->Re(rhos, U0, td.Uc(), this->d_, mus);


    // 3. Compute heat- and momentum transfers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Heat transfer
    // ~~~~~~~~~~~~~

    // Calculate new particle temperature
    this->T_ =
        this->calcHeatTransfer
        (
            cloud,
            td,
            dt,
            Res,
            Prs,
            kappas,
            NCpW,
            Sh,
            dhsTrans,
            Sph
        );


    this->Cp_ = CpEff(cloud, td, pc, this->T_, idG, idL, idS);

    // =========================================================================
    // --- METAL PHASE TRANSITION & MELTING STATE MACHINE (Fig. 4 Validation) --
    // =========================================================================
    {
        const scalar TmeltAl = 933.0;
        const scalar hfusAl = 3.97e5;     // Latent heat of fusion Al [J/kg]
        const scalar TmeltOxide = 2327.0;
        const scalar hfusOxide = 1.07e6;  // Latent heat of fusion Al2O3 [J/kg]

        // 1. Resolve species local indices
        label idAl_s = -1;
        label idAl2O3_s = -1;
        if (idS >= 0)
        {
            idAl_s = composition.localId(idS, "Al", true);
            if (idAl_s < 0) idAl_s = composition.localId(idS, "C", true);
            if (idAl_s < 0 && YSolid_.size() > 0) idAl_s = 0;

            idAl2O3_s = composition.localId(idS, "Al2O3", true);
            if (idAl2O3_s < 0) idAl2O3_s = composition.localId(idS, "ash", true);
            if (idAl2O3_s < 0 && YSolid_.size() > 1) idAl2O3_s = 1;
        }

        label idAl_l = -1;
        if (idL >= 0)
        {
            idAl_l = composition.localId(idL, "Al_liquid", true);
            if (idAl_l < 0) idAl_l = composition.localId(idL, "Al", true);
            if (idAl_l < 0 && YLiquid_.size() > 0) idAl_l = 0;
        }

        // 2. Component masses in parcel
        scalar m_Al_s = 0.0;
        scalar m_ash = 0.0;
        if (idS >= 0)
        {
            if (idAl_s >= 0 && idAl_s < YSolid_.size())
                m_Al_s = mass1 * this->Y_[idS] * YSolid_[idAl_s];
            if (idAl2O3_s >= 0 && idAl2O3_s < YSolid_.size())
                m_ash = mass1 * this->Y_[idS] * YSolid_[idAl2O3_s];
        }

        scalar m_Al_l = 0.0;
        if (idL >= 0 && idAl_l >= 0 && idAl_l < YLiquid_.size())
        {
            m_Al_l = mass1 * this->Y_[idL] * YLiquid_[idAl_l];
        }

        // 3. Aluminum Melting Plateau at 933 K
        if (m_Al_s > 1e-25 && this->T_ >= TmeltAl && idAl_l >= 0)
        {
            scalar E_excess = mass1 * this->Cp_ * (this->T_ - TmeltAl);
            if (E_excess > 0.0)
            {
                scalar dMelt = E_excess / hfusAl;
                if (dMelt < m_Al_s)
                {
                    m_Al_s -= dMelt;
                    m_Al_l += dMelt;
                    this->T_ = TmeltAl;
                }
                else
                {
                    scalar E_rem = E_excess - m_Al_s * hfusAl;
                    m_Al_l += m_Al_s;
                    m_Al_s = 0.0;
                    this->T_ = TmeltAl + E_rem / max(mass1 * this->Cp_, rootVSmall);
                }
            }
        }

        // 4. Alumina Melting Plateau at 2327 K
        // canCombust_ == 1: oxide shell is solid
        // canCombust_ == 2: oxide shell has melted into an oxide cap
        static scalar m_ash_melted = 0.0;
        scalar Y_Oxides = 0.0;
        scalar Y_Oxidel = 0.0;

        if (m_ash > 1e-25)
        {
            if (this->canCombust_ < 2)
            {
                if (this->T_ >= TmeltOxide)
                {
                    scalar E_excess = mass1 * this->Cp_ * (this->T_ - TmeltOxide);
                    scalar m_unmelted = max(m_ash - m_ash_melted, 0.0);
                    scalar dMelt = E_excess / hfusOxide;

                    if (dMelt < m_unmelted)
                    {
                        // Incomplete melting: temperature locked at 2327 K
                        m_ash_melted += dMelt;
                        this->T_ = TmeltOxide;
                        const scalar fMelt = min(m_ash_melted / max(m_ash, rootVSmall), 1.0);
                        Y_Oxides = (1.0 - fMelt) * (m_ash / max(mass1, 1e-25));
                        Y_Oxidel = fMelt * (m_ash / max(mass1, 1e-25));
                    }
                    else
                    {
                        // Solid oxide completely melted into cap!
                        this->canCombust_ = 2;
                        scalar E_rem = E_excess - m_unmelted * hfusOxide;
                        this->T_ = TmeltOxide + E_rem / max(mass1 * this->Cp_, rootVSmall);
                        m_ash_melted = m_ash;
                        Y_Oxides = 0.0;
                        Y_Oxidel = m_ash / max(mass1, 1e-25);
                    }
                }
                else
                {
                    Y_Oxides = m_ash / max(mass1, 1e-25);
                    Y_Oxidel = 0.0;
                }
            }
            else
            {
                // Already molten oxide cap
                Y_Oxides = 0.0;
                Y_Oxidel = m_ash / max(mass1, 1e-25);
            }
        }

        // Thermal equilibrium limit during intense evaporation (Paper Section 4.1.1):
        // Liquid aluminum cannot exceed its boiling point (2792 K) at 1 atm
        const scalar TboilAl = 2792.0;
        if (m_Al_l > 1e-25 && this->T_ > TboilAl)
        {
            this->T_ = TboilAl;
        }

        // 5. Update phase mass fractions
        scalar mSolidNew = m_Al_s + m_ash;
        scalar mLiquidNew = m_Al_l;
        scalar mGasNew = (idG >= 0) ? (mass1 * this->Y_[idG]) : 0.0;
        scalar mTotalNew = max(mSolidNew + mLiquidNew + mGasNew, rootVSmall);

        if (idG >= 0) this->Y_[idG] = mGasNew / mTotalNew;
        if (idL >= 0) this->Y_[idL] = mLiquidNew / mTotalNew;
        if (idS >= 0) this->Y_[idS] = max(1.0 - this->Y_[idG] - this->Y_[idL], 0.0);

        if (idS >= 0 && mSolidNew > 1e-25)
        {
            if (idAl_s >= 0 && idAl_s < YSolid_.size())
                YSolid_[idAl_s] = m_Al_s / mSolidNew;
            if (idAl2O3_s >= 0 && idAl2O3_s < YSolid_.size())
                YSolid_[idAl2O3_s] = m_ash / mSolidNew;
        }
        if (idL >= 0 && mLiquidNew > 1e-25)
        {
            if (idAl_l >= 0 && idAl_l < YLiquid_.size())
                YLiquid_[idAl_l] = 1.0;
        }

        // 6. Two-material composite density and diameter (Eq. 13 of paper)
        const scalar rhoAl = 2700.0;
        const scalar rhoOxide = 3950.0;
        const scalar mAlTot = m_Al_s + m_Al_l;
        const scalar V_p = (mAlTot / rhoAl) + (m_ash / rhoOxide);
        if (V_p > 1e-25)
        {
            this->d_ = cbrt(V_p * 6.0 / constant::mathematical::pi);
            this->rho_ = mass1 / V_p;
        }

        // 7. Refresh Cp
        this->Cp_ = CpEff(cloud, td, pc, this->T_, idG, idL, idS);

        // 8. Fig. 4 High-Resolution Data Logging
        const scalar time_ms = td.mesh.time().value() * 1000.0;
        const scalar mass0_ref = max(this->mass0_, 1.4137e-12);
        const scalar d0_ref = cbrt(mass0_ref / 2700.0 * 6.0 / constant::mathematical::pi);

        scalar dMassAl_evap = 0.0;
        if (idAl_l >= 0 && idAl_l < dMassPC.size())
        {
            dMassAl_evap = max(dMassPC[idAl_l], scalar(0.0));
        }

        scalar dMassAl_HSR = 0.0;
        if (idAl_l >= 0 && idAl_l < dMassSRLiquid.size())
        {
            dMassAl_HSR += max(dMassSRLiquid[idAl_l], scalar(0.0));
        }
        if (idAl_s >= 0 && idAl_s < dMassSRSolid.size())
        {
            dMassAl_HSR += max(dMassSRSolid[idAl_s], scalar(0.0));
        }

        const scalar htc = cloud.heatTransfer().htc(d0, Res, Prs, kappas, NCpW);
        const scalar Q_inter = htc * constant::mathematical::pi * sqr(d0) * (td.Tc() - T0);
        const scalar Q_HSR = (dMassAl_HSR / dt) * 3.1e7;
        const scalar Q_evap = (dMassAl_evap / dt) * 1.05e7;
        scalar Q_rad = 0.0;
        if (cloud.radiation())
        {
            const tetIndices tetIs = this->currentTetIndices(td.mesh);
            const scalar Gc = td.GInterp().interpolate(this->coordinates(), tetIs);
            const scalar sigma = 5.670374419e-8;
            const scalar epsilon = cloud.constProps().epsilon0();
            Q_rad = constant::mathematical::pi * sqr(d0) * epsilon * (Gc/4.0 - sigma*pow4(T0));
        }

        Info<< "FIG4_DATA: time_ms=" << time_ms
            << " Tp=" << this->T_
            << " mp_ratio=" << (mass1 / mass0_ref)
            << " dp_ratio=" << (this->d_ / d0_ref)
            << " Y_Als=" << (m_Al_s / max(mass1, 1e-25))
            << " Y_All=" << (m_Al_l / max(mass1, 1e-25))
            << " Y_Oxides=" << Y_Oxides
            << " Y_Oxidel=" << Y_Oxidel
            << " mdot_HSR=" << (dMassAl_HSR / dt)
            << " mdot_evap=" << (dMassAl_evap / dt)
            << " Q_inter=" << Q_inter
            << " Q_HSR=" << Q_HSR
            << " Q_evap=" << Q_evap
            << " Q_rad=" << Q_rad
            << endl;
    }

    // Motion
    // ~~~~~~

    // Calculate new particle velocity
    this->U_ =
        this->calcVelocity(cloud, td, dt, Res, mus, mass1, Su, dUTrans, Spu);


    // 4. Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (cloud.solution().coupled())
    {
        // Transfer mass lost to carrier mass, momentum and enthalpy sources
        forAll(YGas_, i)
        {
            scalar dm = np0*dMassGas[i];
            label gid = composition.localToCarrierId(idG, i);
            scalar hs = composition.carrier().hsi(gid, pc, T0);
            cloud.rhoTrans(gid)[this->cell()] += dm;
            cloud.UTransRef()[this->cell()] += dm*U0;
            cloud.hsTransRef()[this->cell()] += dm*hs;
        }
        forAll(YLiquid_, i)
        {
            scalar dm = np0*dMassLiquid[i];
            label gid = composition.localToCarrierId(idL, i);
            scalar hs = composition.carrier().hsi(gid, pc, T0);
            cloud.rhoTrans(gid)[this->cell()] += dm;
            cloud.UTransRef()[this->cell()] += dm*U0;
            cloud.hsTransRef()[this->cell()] += dm*hs;
        }

        // No mapping between solid components and carrier phase
        /*
        forAll(YSolid_, i)
        {
            scalar dm = np0*dMassSolid[i];
            label gid = composition.localToCarrierId(idS, i);
            scalar hs = composition.carrier().hsi(gid, pc, T0);
            cloud.rhoTrans(gid)[this->cell()] += dm;
            cloud.UTransRef()[this->cell()] += dm*U0;
            cloud.hsTransRef()[this->cell()] += dm*hs;
        }
        */

        forAll(dMassSRCarrier, i)
        {
            scalar dm = np0*dMassSRCarrier[i];
            scalar hs = composition.carrier().hsi(i, pc, T0);
            cloud.rhoTrans(i)[this->cell()] += dm;
            cloud.UTransRef()[this->cell()] += dm*U0;
            cloud.hsTransRef()[this->cell()] += dm*hs;
        }

        // Update momentum transfer
        cloud.UTransRef()[this->cell()] += np0*dUTrans;
        cloud.UCoeffRef()[this->cell()] += np0*Spu;

        // Update sensible enthalpy transfer
        cloud.hsTransRef()[this->cell()] += np0*dhsTrans;
        cloud.hsCoeffRef()[this->cell()] += np0*Sph;

        // Update radiation fields
        if (cloud.radiation())
        {
            const scalar ap = this->areaP();
            const scalar T4 = pow4(T0);
            cloud.radAreaP()[this->cell()] += dt*np0*ap;
            cloud.radT4()[this->cell()] += dt*np0*T4;
            cloud.radAreaPT4()[this->cell()] += dt*np0*ap*T4;
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::calcDevolatilisation
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const scalar Ts,
    const scalar d,
    const scalar T,
    const scalar mass,
    const scalar mass0,
    const scalarField& YGasEff,
    const scalarField& YLiquidEff,
    const scalarField& YSolidEff,
    label& canCombust,
    scalarField& dMassDV,
    scalar& Sh,
    scalar& N,
    scalar& NCpW,
    scalarField& Cs
) const
{
    // Check that model is active
    if
    (
        isType
        <
            NoDevolatilisation
            <
                typename TrackCloudType::reactingMultiphaseCloudType
            >
        >
        (
            cloud.devolatilisation()
        )
    )
    {
        if (canCombust != -1)
        {
            canCombust = 1;
        }
        return;
    }

    // Initialise demand-driven constants
    (void)cloud.constProps().TDevol();
    (void)cloud.constProps().LDevol();

    // Check that the parcel temperature is within necessary limits for
    // devolatilisation to occur
    if (T < cloud.constProps().TDevol() || canCombust == -1)
    {
        return;
    }

    typedef typename TrackCloudType::thermoCloudType thermoCloudType;
    const CompositionModel<thermoCloudType>& composition =
        cloud.composition();

    const label idG = composition.idGas();

    const typename TrackCloudType::parcelType& p =
        static_cast<const typename TrackCloudType::parcelType&>(*this);
    typename TrackCloudType::parcelType::trackingData& ttd =
        static_cast<typename TrackCloudType::parcelType::trackingData&>(td);

    // Total mass of volatiles evolved
    cloud.devolatilisation().calculate
    (
        p,
        ttd,
        dt,
        mass0,
        mass,
        T,
        YGasEff,
        YLiquidEff,
        YSolidEff,
        canCombust,
        dMassDV
    );

    scalar dMassTot = sum(dMassDV);

    cloud.devolatilisation().addToDevolatilisationMass
    (
        this->nParticle_*dMassTot
    );

    Sh -= dMassTot*cloud.constProps().LDevol()/dt;

    // Update molar emissions
    if (cloud.heatTransfer().BirdCorrection())
    {
        // Molar average molecular weight of carrier mix
        const scalar Wc = max(small, td.rhoc()*RR*td.Tc()/td.pc());

        // Note: hardcoded gaseous diffusivities for now
        // TODO: add to carrier thermo
        const scalar beta = sqr(cbrt(15.0) + cbrt(15.0));

        forAll(dMassDV, i)
        {
            const label id = composition.localToCarrierId(idG, i);
            const scalar Cp = composition.carrier().Cpi(id, td.pc(), Ts);
            const scalar W = composition.carrier().WiValue(id);
            const scalar Ni = dMassDV[i]/(this->areaS(d)*dt*W);

            // Dab calc'd using API vapour mass diffusivity function
            const scalar Dab =
                3.6059e-3*(pow(1.8*Ts, 1.75))
               *sqrt(1.0/W + 1.0/Wc)
               /(td.pc()*beta);

            N += Ni;
            NCpW += Ni*Cp*W;
            Cs[id] += Ni*d/(2.0*Dab);
        }
    }
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseParcel<ParcelType>::calcSurfaceReactions
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const scalar d,
    const scalar T,
    const scalar mass,
    const label canCombust,
    const scalar N,
    const scalarField& YMix,
    const scalarField& YGas,
    const scalarField& YLiquid,
    const scalarField& YSolid,
    scalarField& dMassSRGas,
    scalarField& dMassSRLiquid,
    scalarField& dMassSRSolid,
    scalarField& dMassSRCarrier,
    scalar& Sh,
    scalar& dhsTrans
) const
{
    // Check that model is active
    if
    (
        isType
        <
            NoSurfaceReaction
            <
                typename TrackCloudType::reactingMultiphaseCloudType
            >
        >
        (
            cloud.surfaceReaction()
        )
    )
    {
        return;
    }

    // Initialise demand-driven constants
    (void)cloud.constProps().hRetentionCoeff();
    (void)cloud.constProps().TMax();

    // Check that model is active
    if (canCombust < 1)
    {
        return;
    }


    // Update surface reactions
    const scalar hReaction = cloud.surfaceReaction().calculate
    (
        dt,
        this->cell(),
        d,
        T,
        td.Tc(),
        td.pc(),
        td.rhoc(),
        mass,
        YGas,
        YLiquid,
        YSolid,
        YMix,
        N,
        dMassSRGas,
        dMassSRLiquid,
        dMassSRSolid,
        dMassSRCarrier
    );

    cloud.surfaceReaction().addToSurfaceReactionMass
    (
        this->nParticle_
       *(sum(dMassSRGas) + sum(dMassSRLiquid) + sum(dMassSRSolid))
    );

    const scalar xsi = min(T/cloud.constProps().TMax(), 1.0);
    const scalar coeff =
        (1.0 - xsi*xsi)*cloud.constProps().hRetentionCoeff();

    Sh += coeff*hReaction/dt;

    dhsTrans += (1.0 - coeff)*hReaction;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingMultiphaseParcel<ParcelType>::ReactingMultiphaseParcel
(
    const ReactingMultiphaseParcel<ParcelType>& p
)
:
    ParcelType(p),
    mass0_(p.mass0_),
    YGas_(p.YGas_),
    YLiquid_(p.YLiquid_),
    YSolid_(p.YSolid_),
    canCombust_(p.canCombust_)
{}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "ReactingMultiphaseParcelIO.C"

// ************************************************************************* //
