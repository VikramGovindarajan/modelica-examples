within Interfaces;

partial model PartialLumpedVolume
"Lumped volume with mass and energy balance"
import Types;
import Types.Dynamics;
import Modelica.Media.Interfaces.Choices.IndependentVariables;
import SI=Modelica.Units.SI;

  outer System system "System properties";
  replaceable package Medium =
    Modelica.Media.Interfaces.PartialMedium "Medium in the component";

  // Inputs provided to the volume model
  input SI.Volume fluidVolume "Volume";

  // Assumptions
  parameter Types.Dynamics energyDynamics=system.energyDynamics
  "Formulation of energy balance";

  parameter Types.Dynamics massDynamics=system.massDynamics
  "Formulation of mass balance";

  final parameter Types.Dynamics substanceDynamics=massDynamics
  "Formulation of substance balance";

  final parameter Types.Dynamics traceDynamics=massDynamics
  "Formulation of trace substance balance";

  // Initialization
  parameter Medium.AbsolutePressure p_start = system.p_start
  "Start value of pressure";

  parameter Boolean use_T_start = true
  "= true, use T_start, otherwise h_start";

  parameter Medium.Temperature T_start=
    if use_T_start then system.T_start else Medium.temperature_phX(p_start,h_start,X_start)
  "Start value of temperature";

  parameter Medium.SpecificEnthalpy h_start=
    if use_T_start then Medium.specificEnthalpy_pTX(p_start, T_start, X_start) else Medium.h_default
  "Start value of specific enthalpy";

  parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default
  "Start value of mass fractions m_i/m";
  parameter Medium.ExtraProperty C_start[Medium.nC](
       quantity=Medium.extraPropertiesNames) = Medium.C_default
  "Start value of trace substances";

  Medium.BaseProperties medium(
    preferredMediumStates = (if energyDynamics == Dynamics.SteadyState and
                                massDynamics   == Dynamics.SteadyState then false else true),
    p(start=p_start),
    h(start=h_start),
    T(start=T_start),
    Xi(start=X_start[1:Medium.nXi]));
  SI.Energy U "Internal energy of fluid";
  SI.Mass m "Mass of fluid";
  SI.Mass[Medium.nXi] mXi "Masses of independent components in the fluid";
  SI.Mass[Medium.nC] mC "Masses of trace substances in the fluid";
  // C need to be added here because unlike for Xi, which has medium.Xi,
  // there is no variable medium.C
  Medium.ExtraProperty C[Medium.nC] "Trace substance mixture content";

  // variables that need to be defined by an extending class
  SI.MassFlowRate mb_flow "Mass flows across boundaries";
  SI.MassFlowRate[Medium.nXi] mbXi_flow
  "Substance mass flows across boundaries";
  Medium.ExtraPropertyFlowRate[Medium.nC] mbC_flow
  "Trace substance mass flows across boundaries";
  SI.EnthalpyFlowRate Hb_flow
  "Enthalpy flow across boundaries or energy source/sink";
  SI.HeatFlowRate Qb_flow
  "Heat flow across boundaries or energy source/sink";
  SI.Power Wb_flow "Work flow across boundaries or source term";

protected

  parameter Boolean initialize_p = not Medium.singleState
  "= true to set up initial equations for pressure";

  Real[Medium.nC] mC_scaled(min=fill(Modelica.Constants.eps, Medium.nC))
  "Scaled masses of trace substances in the fluid";

equation
  assert(not (energyDynamics<>Dynamics.SteadyState and massDynamics==Dynamics.SteadyState) or Medium.singleState,
         "Bad combination of dynamics options and Medium not conserving mass if fluidVolume is fixed.");

  // Total quantities
  m = fluidVolume*medium.d;
  mXi = m*medium.Xi;
  U = m*medium.u;
  mC = m*C;

  // Energy and mass balances
  if energyDynamics == Dynamics.SteadyState then
    0 = Hb_flow + Qb_flow + Wb_flow;
  else
    der(U) = Hb_flow + Qb_flow + Wb_flow;
  end if;

  if massDynamics == Dynamics.SteadyState then
    0 = mb_flow;
  else
    der(m) = mb_flow;
  end if;

  if substanceDynamics == Dynamics.SteadyState then
    zeros(Medium.nXi) = mbXi_flow;
  else
    der(mXi) = mbXi_flow;
  end if;

  if traceDynamics == Dynamics.SteadyState then
    zeros(Medium.nC)  = mbC_flow;
  else
    der(mC_scaled) = mbC_flow./Medium.C_nominal;
  end if;
    mC = mC_scaled.*Medium.C_nominal;

initial equation
  // initialization of balances
  if energyDynamics == Dynamics.FixedInitial then
    if Medium.ThermoStates == IndependentVariables.ph or
       Medium.ThermoStates == IndependentVariables.phX then
       medium.h = h_start;
    else
       medium.T = T_start;
    end if;
  elseif energyDynamics == Dynamics.SteadyStateInitial then

    if Medium.ThermoStates == IndependentVariables.ph or
       Medium.ThermoStates == IndependentVariables.phX then
       der(medium.h) = 0;
    else
       der(medium.T) = 0;
    end if;
  end if;

  if massDynamics == Dynamics.FixedInitial then
    if initialize_p then
      medium.p = p_start;
    end if;
  elseif massDynamics == Dynamics.SteadyStateInitial then
    if initialize_p then
      der(medium.p) = 0;
    end if;
  end if;

  if substanceDynamics == Dynamics.FixedInitial then
    medium.Xi = X_start[1:Medium.nXi];
  elseif substanceDynamics == Dynamics.SteadyStateInitial then
    der(medium.Xi) = zeros(Medium.nXi);
  end if;

  if traceDynamics == Dynamics.FixedInitial then
    mC_scaled = m*C_start[1:Medium.nC]./Medium.C_nominal;
  elseif traceDynamics == Dynamics.SteadyStateInitial then
    der(mC_scaled) = zeros(Medium.nC);
  end if;

  
end PartialLumpedVolume;
