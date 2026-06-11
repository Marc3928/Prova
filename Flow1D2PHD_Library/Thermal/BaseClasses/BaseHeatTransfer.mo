within Flow1D2PHD_Library.HeatTransfer.BaseClasses;

partial model BaseHeatTransfer
    replaceable package Medium = Modelica.Media.Interfaces.PartialTwoPhaseMedium "Medium model";
    parameter Modelica.Units.SI.Length Dhyd;
    parameter Modelica.Units.SI.Length L;
    input Modelica.Units.SI.MassFlowRate w;
    input Medium.ThermodynamicState fluidState;
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port annotation(
    Placement(transformation(origin = {-24, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(extent = {{-10, -10}, {10, 10}})));
    output Modelica.Units.SI.Power Q_flow;
end BaseHeatTransfer;
