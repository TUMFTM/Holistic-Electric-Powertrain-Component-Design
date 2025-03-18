function InterpolateMotorSpecs(dataTable, vehicleWeight)
    % INTERPOLATEMOTORSPECS - Estimates max power and max torque based on vehicle weight.
    %
    % Outputs:
    %   - Console feedback on interpolated motor specifications.

    % Predict max power based on vehicle weight
    predictedPower = generalizedInterpolation(dataTable, 'WeightBeforeTeardown_kg_', 'MaxPower_kW_', vehicleWeight);

    % Predict max torque based on vehicle weight
    predictedTorque = generalizedInterpolation(dataTable, 'WeightBeforeTeardown_kg_', 'MaxTorque_Nm_', vehicleWeight);

    % Display estimated motor specs
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
    disp('🔍  **Interpolated Motor Specifications Based on Vehicle Weight**');
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
    fprintf('🚗 **Vehicle Weight:** %.2f kg\n', vehicleWeight);
    fprintf('⚡ **Estimated Max Power (Interpolated from Database):** %.2f kW\n', predictedPower);
    fprintf('🔩 **Estimated Max Torque (Interpolated from Database):** %.2f Nm\n', predictedTorque);
    disp('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
end