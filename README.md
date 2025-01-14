# Holistic Electric Powertrain Component Design

## Overview
This simulation tool provides a comprehensive approach for designing electric powertrain components in battery electric vehicles (BEVs) during early development phases. It enables manufacturers to establish initial powertrain component designs that meet specific vehicle concept requirements while minimizing interdependencies between development divisions. The tool encompasses the complete powertrain, including high-voltage battery, power electronics, electric machine, and transmission, serving as a foundation for further development.

The model's simulation results and performance characteristics have been validated against reference vehicles through dynamometer testing and teardown analysis. It is specifically designed to support optimization approaches focusing on energy consumption, which is crucial for BEV design. Through its modular architecture, components can be independently developed and evaluated, helping to avoid time and cost-intensive iterations in the development process.

## Key Features
* Modular powertrain component simulation
* Multiple modeling methods per component
* Standardized driving cycle simulation (WLTC, NEDC, etc.)
* Acceleration cycle simulation
* Energy consumption analysis

## Developer
The main developers of this tool are Jan Koloch and Nico Rosenberger (both from the Institute for Automotive Technology, Technical University of Munich).

The below stated reference is the main documentation for the tool documented in this repository.

There are several other contributors who worked on different modules of the tool. Here follows an overview of the contributors:

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
The simulation is controlled through the main script `Sim_v2.m`.

## Sources
The tool is documented in the following scientific publication:

Rosenberger, N.; Deininger, S.; Koloch, J.; Lienkamp, M. Holistic Electric Powertrain Component Design for Battery Electric Vehicles in an early Development Phase. World Electr. Veh. J. 2024, 1, 0. https://doi.org/
