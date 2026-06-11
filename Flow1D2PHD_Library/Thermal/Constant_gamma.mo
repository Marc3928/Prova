within Flow1D2PHD_Library.HeatTransfer;

model Constant_gamma
  extends Flow1D2PHD_Library.HeatTransfer.BaseClasses.BaseHeatTransfer;
  parameter Modelica.Units.SI.CoefficientOfHeatTransfer gamma;
equation
  Q_flow = gamma * L * Dhyd * (port.T - Medium.temperature(fluidState)); //potenza termica scambiata 
  port.Q_flow + Q_flow = 0; //non c'è accumulo di potenza 

end Constant_gamma;
