within Vessels.BaseClasses.HeatTransfer;

model IdealHeatTransfer
    "IdealHeatTransfer: Ideal heat transfer without thermal resistance"
  extends PartialVesselHeatTransfer;

equation
  Ts = heatPorts.T;

end IdealHeatTransfer;
