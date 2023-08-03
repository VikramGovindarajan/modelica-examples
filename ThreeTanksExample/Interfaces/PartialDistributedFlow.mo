within Interfaces;

partial model PartialDistributedFlow
"Base class for a distributed momentum balance"

  outer System system "System properties";

  replaceable package Medium =
    Modelica.Media.Interfaces.PartialMedium "Medium in the component";
  
  import SI=Modelica.Units.SI;
  
  parameter Boolean allowFlowReversal = system.allowFlowReversal
"= true to allow flow reversal, false restricts to design direction (m_flows >= zeros(m))";

  // Discretization
  parameter Integer m=1 "Number of flow segments";

  // Inputs provided to the flow model
  input SI.Length[m] pathLengths "Lengths along flow path";

  // Variables defined by momentum model
  Medium.MassFlowRate[m] m_flows(
     each min=if allowFlowReversal then -Modelica.Constants.inf else 0,
     each start = m_flow_start,
     each stateSelect = if momentumDynamics == Types.Dynamics.SteadyState then StateSelect.default else
                               StateSelect.prefer)
"Mass flow rates between states";

  // Parameters
  parameter Types.Dynamics momentumDynamics=system.momentumDynamics
"Formulation of momentum balance";

  parameter Medium.MassFlowRate m_flow_start=system.m_flow_start
"Start value of mass flow rates";

  // Total quantities
  SI.Momentum[m] Is "Momenta of flow segments";

  // Source terms and forces to be defined by an extending model (zero if not used)
  SI.Force[m] Ib_flows "Flow of momentum across boundaries";
  SI.Force[m] Fs_p "Pressure forces";
  SI.Force[m] Fs_fg "Friction and gravity forces";

equation
  // Total quantities
  Is = {m_flows[i]*pathLengths[i] for i in 1:m};

  // Momentum balances
  if momentumDynamics == Types.Dynamics.SteadyState then
    zeros(m) = Ib_flows - Fs_p - Fs_fg;
  else
    der(Is) = Ib_flows - Fs_p - Fs_fg;
  end if;

initial equation
  if momentumDynamics == Types.Dynamics.FixedInitial then
    m_flows = fill(m_flow_start, m);
  elseif momentumDynamics == Types.Dynamics.SteadyStateInitial then
    der(m_flows) = zeros(m);
  end if;

end PartialDistributedFlow;
