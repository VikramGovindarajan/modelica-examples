within;
model System
  "System properties and default values (ambient, flow direction, initialization)"

  package Medium = Modelica.Media.Interfaces.PartialMedium
    "Medium model for default start values";
  
  import SI=Modelica.Units.SI;
  
  parameter SI.AbsolutePressure p_ambient=101325
    "Default ambient pressure";
  parameter SI.Temperature T_ambient=293.15
    "Default ambient temperature";
  parameter SI.Acceleration g=Modelica.Constants.g_n
    "Constant gravity acceleration";

  // Assumptions
  parameter Boolean allowFlowReversal = true
    "= false to restrict to design flow direction (port_a -> port_b)";
  parameter Types.Dynamics energyDynamics=
    Types.Dynamics.DynamicFreeInitial
    "Default formulation of energy balances";
  parameter Types.Dynamics massDynamics=
    energyDynamics "Default formulation of mass balances";
  final parameter Types.Dynamics substanceDynamics=
    massDynamics "Default formulation of substance balances";
  final parameter Types.Dynamics traceDynamics=
    massDynamics "Default formulation of trace substance balances";
  parameter Types.Dynamics momentumDynamics=
    Types.Dynamics.SteadyState
    "Default formulation of momentum balances, if options available";

  // Initialization
  parameter SI.MassFlowRate m_flow_start = 0
    "Default start value for mass flow rates";
  parameter SI.AbsolutePressure p_start = p_ambient
    "Default start value for pressures";
  parameter SI.Temperature T_start = T_ambient
    "Default start value for temperatures";
  // Advanced
  parameter Boolean use_eps_Re = false
    "= true to determine turbulent region automatically using Reynolds number";
  parameter SI.MassFlowRate m_flow_nominal = if use_eps_Re then 1 else 1e2*m_flow_small
    "Default nominal mass flow rate";
  parameter Real eps_m_flow(min=0) = 1e-4
    "Regularization of zero flow for |m_flow| < eps_m_flow*m_flow_nominal";
  parameter SI.AbsolutePressure dp_small(min=0) = 1
    "Default small pressure drop for regularization of laminar and zero flow";
  parameter SI.MassFlowRate m_flow_small(min=0) = 1e-2
    "Default small mass flow rate for regularization of laminar and zero flow";
initial equation
  //assert(use_eps_Re, "*** Using classic system.m_flow_small and system.dp_small."
  //       + " They do not distinguish between laminar flow and regularization of zero flow."
  //       + " Absolute small values are error prone for models with local nominal values."
  //       + " Moreover dp_small can generally be obtained automatically."
  //       + " Please update the model to new system.use_eps_Re = true  (see system, Advanced tab). ***",
  //       level=AssertionLevel.warning);
  
end System;
