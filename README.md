# Holistic Electric Powertrain Component Design

## Overview
This simulation tool offers a modular and flexible approach for designing electric powertrain components in battery electric vehicles (BEVs), with an emphasis on performance and energy consumption. This version introduces improved usability, expandability, and simulation accuracy, providing manufacturers with a powerful tool for optimizing powertrain components during early development stages.

The tool simulates key powertrain components, including the high-voltage battery, power electronics, electric motor, and transmission. With user guidance for independent configuration and validation of each component, it reduces iterations and streamlines the development process. The simulation tool also supports optimization calculations for both energy consumption and performance.

## Key Features
* User-friendly framework
* Modular and independent component configuration for powertrain
* Performance and consumption analysis
* Component validation stations
* Parameter dynamics tool to understand parameter interdependencies
* User guide for design suggestions

## Developer
The main developers of this tool are Jan Koloch and Nico Rosenberger (both from the Institute for Automotive Technology, Technical University of Munich).

The below stated reference is the main documentation for the tool documented in this repository.

There are several other contributors who worked on different modules of the tool. Here follows an overview of the contributors:

* Shehata, Poula (Master's Thesis at the Technical University of Munich)
* Sanftl, Stephan (Semester Thesis at the Technical University of Munich)
* Elbadawi, Abdelrahman (Semester Thesis at the Technical University of Munich)

## Requirements
* MATLAB/Simulink (Version 2023a or newer)
* Simscape Battery (for battery modeling)

## Installation
1. Install MATLAB/Simulink and required toolboxes
2. Clone this repository
3. Place the tool folder in your preferred working directory

## Usage
The simulation is controlled through the main script `PowertrainFramework.mlx`.

## Sources
The tool is documented in the following scientific publication:
Rosenberger, N.; Deininger, S.; Koloch, J.; Lienkamp, M. Holistic Electric Powertrain Component Design for Battery Electric Vehicles in an early Development Phase. World Electr. Veh. J. 2024, 1, 0. https://doi.org/10.3390/wevj16020061

