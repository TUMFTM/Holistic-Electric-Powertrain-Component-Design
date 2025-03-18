function cost = calcBatteryCost(battery_GD)

   cost_per_kWh = 157;
   cost = battery_GD.SysInfo.E_sys * cost_per_kWh ;

end