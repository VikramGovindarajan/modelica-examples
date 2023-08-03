within Pipes.BaseClasses.FlowModels;

partial model PartialGenericPipeFlow
"GenericPipeFlow: Pipe flow pressure loss and gravity with replaceable WallFriction package"

  parameter Boolean from_dp = momentumDynamics >= Types.Dynamics.SteadyStateInitial
"= true, use m_flow = f(dp), otherwise dp = f(m_flow)";

  extends
Pipes.BaseClasses.FlowModels.PartialStaggeredFlowModel(
 final Re_turbulent=4000);

  replaceable package WallFriction =
    Pipes.BaseClasses.WallFriction.Detailed
      constrainedby
Pipes.BaseClasses.WallFriction.PartialWallFriction
"Wall friction model";

  import SI=Modelica.Units.SI;
  
  input SI.Length[n-1] pathLengths_internal
"pathLengths used internally; to be defined by extending class";
  input SI.ReynoldsNumber[n-1] Res_turbulent_internal = Re_turbulent*ones(n-1)
"Re_turbulent used internally; to be defined by extending class";

  // Parameters
  parameter SI.AbsolutePressure dp_nominal
"Nominal pressure loss (only for nominal models)";
  parameter SI.MassFlowRate m_flow_nominal "Nominal mass flow rate";
  parameter SI.MassFlowRate m_flow_small = if system.use_eps_Re then system.eps_m_flow*m_flow_nominal else system.m_flow_small
"Within regularization if |m_flows| < m_flow_small (may be wider for large discontinuities in static head)";

protected
  parameter SI.AbsolutePressure dp_small(start = 1, fixed = false)
"Within regularization if |dp| < dp_small (may be wider for large discontinuities in static head)";
  final parameter Boolean constantPressureLossCoefficient=
     use_rho_nominal and (use_mu_nominal or not WallFriction.use_mu)
"= true, if the pressure loss does not depend on fluid states";
  final parameter Boolean continuousFlowReversal=
     (not useUpstreamScheme)
     or constantPressureLossCoefficient
     or not allowFlowReversal
"= true, if the pressure loss is continuous around zero flow";

  SI.Length[n-1] diameters = 0.5*(dimensions[1:n-1] + dimensions[2:n])
"Mean diameters between segments";
  SI.AbsolutePressure dp_fric_nominal=
    sum(WallFriction.pressureLoss_m_flow(
                   m_flow_nominal/nParallel,
                   rho_nominal,
                   rho_nominal,
                   mu_nominal,
                   mu_nominal,
                   pathLengths_internal,
                   diameters,
                   (crossAreas[1:n-1]+crossAreas[2:n])/2,
                   (roughnesses[1:n-1]+roughnesses[2:n])/2,
                   m_flow_small/nParallel,
                   Res_turbulent_internal))
"Pressure loss for nominal conditions";

initial equation
  // initialize dp_small from flow model
  if system.use_eps_Re then
    dp_small = dp_fric_nominal/m_flow_nominal*m_flow_small;
  else
    dp_small = system.dp_small;
  end if;

equation
  for i in 1:n-1 loop
    assert(m_flows[i] > -m_flow_small or allowFlowReversal, "Reversing flow occurs even though allowFlowReversal is false");
  end for;

  if continuousFlowReversal then
    // simple regularization
    if from_dp and not WallFriction.dp_is_zero then
      m_flows = homotopy(
        actual=  WallFriction.massFlowRate_dp(
                   dps_fg - {g*dheights[i]*rhos_act[i] for i in 1:n-1},
                   rhos_act,
                   rhos_act,
                   mus_act,
                   mus_act,
                   pathLengths_internal,
                   diameters,
                   (crossAreas[1:n-1]+crossAreas[2:n])/2,
                   (roughnesses[1:n-1]+roughnesses[2:n])/2,
                   dp_small/(n-1),
                   Res_turbulent_internal)*nParallel,
        simplified=  m_flow_nominal/dp_nominal*(dps_fg - g*dheights*rho_nominal));
    else
      dps_fg = homotopy(
        actual=  WallFriction.pressureLoss_m_flow(
                   m_flows/nParallel,
                   rhos_act,
                   rhos_act,
                   mus_act,
                   mus_act,
                   pathLengths_internal,
                   diameters,
                   (crossAreas[1:n-1]+crossAreas[2:n])/2,
                   (roughnesses[1:n-1]+roughnesses[2:n])/2,
                   m_flow_small/nParallel,
                   Res_turbulent_internal) + {g*dheights[i]*rhos_act[i] for i in 1:n-1},
        simplified=  dp_nominal/m_flow_nominal*m_flows + g*dheights*rho_nominal);
    end if;
  else
    // regularization for discontinuous flow reversal and static head
    if from_dp and not WallFriction.dp_is_zero then
      m_flows = homotopy(
        actual=  WallFriction.massFlowRate_dp_staticHead(
                   dps_fg,
                   rhos[1:n-1],
                   rhos[2:n],
                   mus[1:n-1],
                   mus[2:n],
                   pathLengths_internal,
                   diameters,
                   g*dheights,
                   (crossAreas[1:n-1]+crossAreas[2:n])/2,
                   (roughnesses[1:n-1]+roughnesses[2:n])/2,
                   dp_small/(n-1),
                   Res_turbulent_internal)*nParallel,
        simplified=  m_flow_nominal/dp_nominal*(dps_fg - g*dheights*rho_nominal));
    else
      dps_fg = homotopy(
        actual=  WallFriction.pressureLoss_m_flow_staticHead(
                   m_flows/nParallel,
                   rhos[1:n-1],
                   rhos[2:n],
                   mus[1:n-1],
                   mus[2:n],
                   pathLengths_internal,
                   diameters,
                   g*dheights,
                   (crossAreas[1:n-1]+crossAreas[2:n])/2,
                   (roughnesses[1:n-1]+roughnesses[2:n])/2,
                   m_flow_small/nParallel,
                   Res_turbulent_internal),
        simplified=  dp_nominal/m_flow_nominal*m_flows + g*dheights*rho_nominal);
    end if;
  end if;

end PartialGenericPipeFlow;
