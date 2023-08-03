within Pipes.BaseClasses.FlowModels;

partial model PartialStaggeredFlowModel
"Base class for momentum balances in flow models"

  //
  // Internal interface
  // (not exposed to GUI; needs to be hard coded when using this model
  //
  replaceable package Medium =
    Modelica.Media.Interfaces.PartialMedium "Medium in the component";
  
  import SI=Modelica.Units.SI;

  parameter Integer n=2 "Number of discrete flow volumes";

  // Inputs
  input Medium.ThermodynamicState[n] states
"Thermodynamic states along design flow";
  input SI.Velocity[n] vs
"Mean velocities of fluid flow";

  // Geometry parameters and inputs
  parameter Real nParallel
"Number of identical parallel flow devices";

  input SI.Area[n] crossAreas
"Cross flow areas at segment boundaries";
  input SI.Length[n] dimensions
"Characteristic dimensions for fluid flow (diameters for pipe flow)";
  input Types.Roughness[n] roughnesses
"Average height of surface asperities";

  // Static head
  input SI.Length[n-1] dheights
"Height(states[2:n]) - Height(states[1:n-1])";

  parameter SI.Acceleration g=system.g
"Constant gravity acceleration";

  // Assumptions
  parameter Boolean allowFlowReversal=system.allowFlowReversal
"= true, if flow reversal is enabled, otherwise restrict flow to design direction (states[1] -> states[n+1])";
  parameter Types.Dynamics momentumDynamics=system.momentumDynamics
"Formulation of momentum balance";

  // Initialization
  parameter Medium.MassFlowRate m_flow_start=system.m_flow_start
"Start value of mass flow rates";
  parameter Medium.AbsolutePressure p_a_start
"Start value for p[1] at design inflow";
  parameter Medium.AbsolutePressure p_b_start
"Start value for p[n+1] at design outflow";

  //
  // Implementation of momentum balance
  //
  extends Interfaces.PartialDistributedFlow(
             final m = n-1);

  // Advanced parameters
  parameter Boolean useUpstreamScheme = true
"= false to average upstream and downstream properties across flow segments";

  parameter Boolean use_Ib_flows = momentumDynamics <> Types.Dynamics.SteadyState
"= true to consider differences in flow of momentum through boundaries";

  // Variables
  Medium.Density[n] rhos = if use_rho_nominal then fill(rho_nominal, n) else Medium.density(states);
  Medium.Density[n-1] rhos_act "Actual density per segment";

  Medium.DynamicViscosity[n] mus = if use_mu_nominal then fill(mu_nominal, n) else Medium.dynamicViscosity(states);
  Medium.DynamicViscosity[n-1] mus_act "Actual viscosity per segment";

  // Variables
  SI.Pressure[n-1] dps_fg(each start = (p_a_start - p_b_start)/(n-1))
"Pressure drop between states";

  // Reynolds Number
  parameter SI.ReynoldsNumber Re_turbulent = 4000
"Start of turbulent regime, depending on type of flow device";
  parameter Boolean show_Res = false
"= true, if Reynolds numbers are included for plotting";
  SI.ReynoldsNumber[n] Res=Pipes.BaseClasses.CharacteristicNumbers.ReynoldsNumber(
      vs,
      rhos,
      mus,
      dimensions) if show_Res "Reynolds numbers";
  Medium.MassFlowRate[n-1] m_flows_turbulent=
      {nParallel*(crossAreas[i] + crossAreas[i+1])/(dimensions[i] + dimensions[i+1])*mus_act[i]*Re_turbulent for i in 1:n-1} if
         show_Res "Start of turbulent flow";
protected
  parameter Boolean use_rho_nominal = false
"= true, if rho_nominal is used, otherwise computed from medium";
  parameter SI.Density rho_nominal = Medium.density_pTX(Medium.p_default, Medium.T_default, Medium.X_default)
"Nominal density (e.g., rho_liquidWater = 995, rho_air = 1.2)";

  parameter Boolean use_mu_nominal = false
"= true, if mu_nominal is used, otherwise computed from medium";
  parameter SI.DynamicViscosity mu_nominal = Medium.dynamicViscosity(
                                                 Medium.setState_pTX(
                                                     Medium.p_default, Medium.T_default, Medium.X_default))
"Nominal dynamic viscosity (e.g., mu_liquidWater = 1e-3, mu_air = 1.8e-5)";

equation
  if not allowFlowReversal then
    rhos_act = rhos[1:n-1];
    mus_act = mus[1:n-1];
  elseif not useUpstreamScheme then
    rhos_act = 0.5*(rhos[1:n-1] + rhos[2:n]);
    mus_act = 0.5*(mus[1:n-1] + mus[2:n]);
  else
    for i in 1:n-1 loop
      rhos_act[i] = noEvent(if m_flows[i] > 0 then rhos[i] else rhos[i+1]);
      mus_act[i] = noEvent(if m_flows[i] > 0 then mus[i] else mus[i+1]);
    end for;
  end if;

  if use_Ib_flows then
    Ib_flows = nParallel*{rhos[i]*vs[i]*vs[i]*crossAreas[i] - rhos[i+1]*vs[i+1]*vs[i+1]*crossAreas[i+1] for i in 1:n-1};
    // alternatively use densities rhos_act of actual streams, together with mass flow rates,
    // not conserving momentum if fluid density changes between flow segments:
    //Ib_flows = {((rhos[i]*vs[i])^2*crossAreas[i] - (rhos[i+1]*vs[i+1])^2*crossAreas[i+1])/rhos_act[i] for i in 1:n-1};
  else
    Ib_flows = zeros(n-1);
  end if;

  Fs_p = nParallel*{0.5*(crossAreas[i]+crossAreas[i+1])*(Medium.pressure(states[i+1])-Medium.pressure(states[i])) for i in 1:n-1};

  // Note: the equation is written for dps_fg instead of Fs_fg to help the translator
  dps_fg = {Fs_fg[i]/nParallel*2/(crossAreas[i]+crossAreas[i+1]) for i in 1:n-1};

end PartialStaggeredFlowModel;
