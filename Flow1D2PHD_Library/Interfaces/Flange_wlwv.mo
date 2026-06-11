within Flow1D2PHD_Library.Interfaces;

connector Flange_wlwv
    replaceable package Medium = Modelica.Media.Interfaces.PartialTwoPhaseMedium "Medium model";
    ThermoPower.FluidPh.Flange liq(redeclare package Medium = Medium);
    ThermoPower.FluidPh.Flange vap(redeclare package Medium = Medium);
    Modelica.Units.SI.PerUnit alpha "Void fraction";
end Flange_wlwv;
