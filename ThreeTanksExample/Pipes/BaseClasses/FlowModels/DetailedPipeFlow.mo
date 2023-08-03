within Pipes.BaseClasses.FlowModels;

model DetailedPipeFlow
"DetailedPipeFlow: Detailed characteristic for laminar and turbulent flow"
  extends
Pipes.BaseClasses.FlowModels.PartialGenericPipeFlow(
redeclare package WallFriction =
    Pipes.BaseClasses.WallFriction.Detailed,
pathLengths_internal=pathLengths,
dp_nominal(start=1, fixed=false),
m_flow_nominal=if system.use_eps_Re then system.m_flow_nominal else 1e2*m_flow_small,
Res_turbulent_internal = Re_turbulent*ones(n-1));

initial equation
  // initialize dp_nominal from flow model
  if system.use_eps_Re then
    dp_nominal = dp_fric_nominal + g*sum(dheights)*rho_nominal;
  else
    dp_nominal = 1e3*dp_small;
  end if;

end DetailedPipeFlow;
